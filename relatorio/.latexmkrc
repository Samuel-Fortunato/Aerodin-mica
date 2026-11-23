# Compile using LuaLaTeX
$pdf_mode = 4;

# LuaLaTeX options
$lualatex = 'lualatex -synctex=1 -interaction=nonstopmode -file-line-error -halt-on-error %O %S';

# Keep generated files out of the source directory
$out_dir = 'build';

# Main document
@default_files = ('main.tex');
