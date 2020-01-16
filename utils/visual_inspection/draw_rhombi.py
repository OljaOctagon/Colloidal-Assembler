import cairo
import math

WIDTH = 3 
HEIGHT = 3
PIXEL_SCALE = 100


with cairo.SVGSurface("example.svg", 1000, 1000) as surface:
    ctx = cairo.Context(surface)
    ctx.scale(PIXEL_SCALE, PIXEL_SCALE)

    ctx.rectangle(0, 0, WIDTH, HEIGHT)
    ctx.set_source_rgb(0.8, 0.8, 1)
    ctx.fill()

    ctx.move_to(1, 0.5)
    ctx.line_to(2, 0.5)
    ctx.line_to(2.2, 1.3)
    ctx.line_to(1.5, 1.7)
    ctx.line_to(0.8, 1.3)
    ctx.close_path()

    ctx.set_source_rgb(1, 0.5, 0)
    ctx.fill_preserve()

    ctx.set_source_rgb(1, 1, 0)
    ctx.set_line_width(0.04)
    ctx.stroke()

# Drawing code
# ctx.move_to(1, 1)
# ctx.line_to(2.5, 1.5)

# ctx.set_source_rgb(1, 0, 0)
# ctx.set_line_width(0.06)
# ctx.stroke()
# End of drawing code

#surface.write('polyon.svg')


