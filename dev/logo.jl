# Please run `include("./dev/logo.jl")` on RELP.
# https://convertio.co/ja/svg-ico/

as = 4
lw = 9

R = 125
r = 36

x1 = 162.5 - R/sqrt(2)
x2 = 162.5 + R/sqrt(2)
x3 = 162.5 - R/sqrt(2)
x4 = 162.5 + R/sqrt(2)

y1 = 150.0 + R/sqrt(2)
y2 = 150.0 - R/sqrt(2)
y3 = 150.0 - R/sqrt(2)
y4 = 150.0 + R/sqrt(2)

c1 = "#CB3C33"
c2 = "#4063D8"
c3 = "#389826"
c4 = "#9558B2"

function axy(x1, y1, x2, y2)
  dx = x2 - x1
  dy = y2 - y1
  ds = sqrt(dx^2 + dy^2)
  dx = dx / ds * as * lw
  dy = dy / ds * as * lw
  return "$(x2 - dx),$(y2 - dy)"
end

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

  <!-- Auxiliary -->
  <!-- <rect width="325" height="300" fill="#FFFFFF"/> -->
  <!-- <circle cx="162.5" cy="150.0" r="$(R)" stroke="#24292E" stroke-width="9" stroke-linecap="round" stroke-linejoin="round" fill="none"/> -->

  <!-- Circle -->
  <circle cx="$(x1)" cy="$(y1)" r="$(r)" fill="$(c1)"/>
  <circle cx="$(x2)" cy="$(y2)" r="$(r)" fill="$(c2)"/>
  <!-- <circle cx="$(x3)" cy="$(y3)" r="$(r)" fill="$(c3)"/> -->
  <!-- <circle cx="$(x4)" cy="$(y4)" r="$(r)" fill="$(c4)"/> -->

  <!-- Defs -->
  <defs>
    <marker id="arrow" markerWidth="$(as)" markerHeight="$(as)" refX="0" refY="$(as/2)" orient="auto">
      <polygon points="0,0 0,$(as) $(as),$(as/2)" fill="#24292E"/>
    </marker>
  </defs>

  <!-- Arrow -->
  <path d="M$(x1),$(y1) L$(axy(x1,y1,x2,y2)))" stroke="#24292E" stroke-width="$(lw)" stroke-linecap="round"  marker-end="url(#arrow)"/>
  <!-- <path d="M$(x1),$(y1) L$(axy(x1,y1,x3,y3)))" stroke="#24292E" stroke-width="$(lw)" stroke-linecap="round"  marker-end="url(#arrow)"/> -->
  <!-- <path d="M$(x1),$(y1) L$(axy(x1,y1,x4,y4)))" stroke="#24292E" stroke-width="$(lw)" stroke-linecap="round"  marker-end="url(#arrow)"/> -->

</svg>
"""

# HTML(svg) |> display

path = "./docs/src/assets/logo.svg"
mkpath(dirname(path))
file = open(path, "w")
Base.write(file, svg)
close(file)
