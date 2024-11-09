# Please run `include("./dev/logo.jl")` on RELP.
# https://convertio.co/ja/svg-ico/

r = 27.5

using Printf

svg = """
<?xml version="1.0" encoding="UTF-8"?>
<svg
  version="1.1"
  xmlns="http://www.w3.org/2000/svg"
  xmlns:xlink="http://www.w3.org/1999/xlink"
  width="325pt"
  height="300pt"
  viewBox="0 0 325 300"
>

<!-- <rect width="325" height="300" fill="#FFFFFF"/> -->
<circle cx="162.5" cy="150.0" r="125" stroke="#24292E" stroke-width="9" stroke-linecap="round" stroke-linejoin="round" fill="none"/>
<!-- <circle cx="162.5" cy="150.0" r="$r" fill="#CB3C33"/> -->
<circle cx="$(162.5-(150-r)/sqrt(2))" cy="$(150.0+(150-r)/sqrt(2))" r="$r" fill="#CB3C33"/>
<circle cx="$(162.5+(150-r)/sqrt(2))" cy="$(150.0-(150-r)/sqrt(2))" r="$r" fill="#4063D8"/>

</svg>
"""

# HTML(svg) |> display

path = "./docs/src/assets/logo.svg"
mkpath(dirname(path))
file = open(path, "w")
Base.write(file, svg)
close(file)
