#!/usr/bin/env python3
"""Convert a 2-column airfoil .dat file to an extruded Wavefront OBJ surface.

Supports Selig-style files with an optional name/header line, blank lines and
# comments. Numeric coordinate rows must contain x y (or x y z; z is ignored).
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path


def read_coordinates(path: Path):
    points=[]
    for lineno, raw in enumerate(path.read_text(encoding='utf-8', errors='replace').splitlines(),1):
        line=raw.split('#',1)[0].strip()
        if not line: continue
        fields=line.replace(',', ' ').split()
        try:
            vals=[float(v) for v in fields]
        except ValueError:
            # Typical Selig .dat first line is an airfoil name/header.
            if not points: continue
            raise ValueError(f"Non-numeric coordinate at {path}:{lineno}: {raw}")
        if len(vals)<2:
            raise ValueError(f"Expected x y at {path}:{lineno}: {raw}")
        points.append((vals[0], vals[1]))
    if len(points)<3:
        raise ValueError(f"{path}: at least three coordinate points are required")
    # Avoid a duplicate closing point: the writer closes the perimeter itself.
    if len(points)>3 and math.dist(points[0],points[-1])<1e-12:
        points.pop()
    return points


def rotate(points, angle_deg, origin):
    a=math.radians(angle_deg); c=math.cos(a); s=math.sin(a)
    ox,oy=origin
    return [(ox+(x-ox)*c-(y-oy)*s, oy+(x-ox)*s+(y-oy)*c) for x,y in points]


def write_obj(points, output: Path, span: float, object_name: str, zone_name: str):
    if span<=0: raise ValueError('--span must be positive')
    half=span/2; n=len(points); nf=2*n
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',encoding='utf-8',newline='\n') as f:
        f.write(f"# Generated from airfoil .dat by airfoil.py\no {object_name}\n")
        f.write(f"# points : {2*n}\n# faces : {nf}\n# zones : 1\n# 0 {zone_name} (nFaces: {nf})\n")
        for x,y in points:
            f.write(f"v {x:.12g} {y:.12g} {-half:.12g}\n")
            f.write(f"v {x:.12g} {y:.12g} {half:.12g}\n")
        f.write(f"g {zone_name}\n")
        for i in range(n):
            j=(i+1)%n
            i0,i1,j0,j1=2*i+1,2*i+2,2*j+1,2*j+2
            f.write(f"f {i0} {i1} {j1}\nf {j1} {j0} {i0}\n")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input','-i',required=True,type=Path,help='Airfoil .dat file')
    p.add_argument('--output','-o',required=True,type=Path,help='Output OBJ file')
    p.add_argument('--angle-deg',type=float,default=0.0,help='CCW rotation in x-y plane')
    p.add_argument('--span',type=float,default=0.2,help='Extrusion thickness in z')
    p.add_argument('--origin',nargs=2,type=float,default=(0.0,0.0),metavar=('X','Y'))
    p.add_argument('--object-name',default='wing')
    p.add_argument('--zone-name',default='WALL10')
    a=p.parse_args()
    pts=rotate(read_coordinates(a.input),a.angle_deg,a.origin)
    write_obj(pts,a.output,a.span,a.object_name,a.zone_name)
    print(f"Generated {a.output}: {len(pts)} perimeter points, {2*len(pts)} vertices, {2*len(pts)} triangular faces, span={a.span:g}, angle={a.angle_deg:g} deg")

if __name__=='__main__': main()
