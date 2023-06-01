module Interactions
export Sliders

using HypertextLiteral
using PlutoUI

render_slider(Child, name, range) =
  @htl("""<tr><td>$name</td><td>$(Child(Slider(range)))</td></tr>""")

function Sliders(title, vars)
  PlutoUI.combine() do Child
    @htl("""
<h3>$title</h3>
<table>
<thead>
<tr><td>Variable</td><td>Value</td></tr>
</thead>
<tbody>
$((render_slider(Child, name, range) for (name, range) in vars))
</tbody>
</table>
    """)
  end
end

end
