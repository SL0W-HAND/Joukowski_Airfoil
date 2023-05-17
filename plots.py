from manim import *
import numpy as np


eta0 = complex(0.2, 0.15)
u0 = 100
gamma = 4*np.pi*u0*(eta0.imag)
a = 1.2
R = 1


def circle(theta):
    return np.exp(complex(0, np.pi*theta))


def shiftfunction(function):
    return function + eta0


def joukouski(z):
    return z + ((a**2))/z


def inverseShiftfunction(eta):
    return eta - eta0


def inverseJoukouski1(eta):
    return 0.5*(eta - np.sqrt((eta**2)-4*(a**2)))


def inverseJoukouski2(eta):
    return 0.5*(eta + np.sqrt((eta**2)-4*(a**2)))

    # positive value


class JoukouskiTransform(Scene):
    def construct(self):
        self.camera.background_color = "#FFF"
        # create two Cartesian planes
        plane1 = NumberPlane(
            x_range=(-1.1, 1.1),
            y_range=(-1.1, 1.1),
            background_line_style={
                "stroke_color": WHITE,
                "stroke_width": 0,
                "stroke_opacity": 0.6
            },
            axis_config={"color": BLACK}
        )
        plane2 = NumberPlane(
            x_range=(-1.1, (1.1+eta0.real)),
            y_range=(-1, 1.3, (1.1+eta0.imag)),
            background_line_style={
                "stroke_color": WHITE,
                "stroke_width": 0,
                "stroke_opacity": 0.6
            },
            axis_config={"color": BLACK}
        )
        plane3 = NumberPlane(
            x_range=(-3, 3),
            y_range=(-0.5, 2),
            background_line_style={
                "stroke_color": WHITE,
                "stroke_width": 0,
                "stroke_opacity": 0.6
            },
            axis_config={"color": BLACK}
        )

        # position the axes
        plane1.to_corner(UL).shift(DOWN*0.3)
        plane3.to_corner(UR).shift(DOWN*0.3)
        plane2.to_edge(DOWN)
        axes = VGroup(plane1, plane2, plane3)

        self.add(axes)

        dot = Dot(color="#000").move_to(
            plane2.coords_to_point(eta0.real, eta0.imag))
        dotLabel = Tex(r"$\eta_0$", font_size=30,
                       color="#000").next_to(dot, 0.5 * UR)
        self.add(dotLabel, dot)
        # create the arrow
        arrow = Arrow(
            start=plane1.get_bottom(),
            end=plane2.get_left(),
            buff=SMALL_BUFF,
            color="#000",
        )

        label1 = Tex(r"$\eta=z + \eta_0$", font_size=40,
                     color="#000").next_to(plane2, LEFT*8.5 + UP)

        title1 = Tex(r"$Plano \quad z$", font_size=40,
                     color="#000").next_to(plane1, UP)

        title2 = Tex(r"$Plano \quad \eta$", font_size=40,
                     color="#000").next_to(plane2, UP)

        title3 = Tex(r"$Plano \quad \zeta$", font_size=40,
                     color="#000").next_to(plane3, UP)
        # add the arrow to the scene
        titles = VGroup(title1, title2, title3)
        self.add(titles)

        arrow2 = Arrow(
            start=plane1.get_right(),
            end=plane3.get_left(),
            buff=SMALL_BUFF,
            color="#000"
        )
        label2 = Tex(
            r"$\zeta = z + \eta_0  + \frac{a^2}{z + \eta_0}$", font_size=40, color="#000").next_to(arrow2, UP)

        arrow3 = Arrow(
            start=plane2.get_right(),
            end=plane3.get_bottom(),
            buff=SMALL_BUFF,
            color="#000"
        )
        label3 = Tex(r"$\zeta = \eta+\frac{a^2}{\eta}$",
                     font_size=40, color="#000").next_to(plane2, RIGHT*7+UP)

        self.add(arrow, label1, arrow2, label2, arrow3, label3)

        curve = plane1.plot_parametric_curve(
            lambda t: np.array(
                [
                    circle(t).real,
                    circle(t).imag,
                    0,
                ]
            ),
            t_range=[0, 2 * np.pi],
            color="#000",
        )

        curve2 = plane2.plot_parametric_curve(
            lambda t: np.array(
                [
                    shiftfunction(circle(t)).real,
                    shiftfunction(circle(t)).imag,
                    0,
                ]
            ),
            t_range=[0, 2 * np.pi],
            color="#000",
        )

        curve3 = plane3.plot_parametric_curve(
            lambda t: np.array(
                [
                    joukouski(shiftfunction(circle(t))).real,
                    joukouski(shiftfunction(circle(t))).imag,
                    0,
                ]
            ),
            t_range=[0, 2 * np.pi],
            color="#000",
        )

        self.add(curve, curve2, curve3)


class Airfoil(Scene):
    def construct(self):
        self.camera.background_color = "#FFF"
        sq = Square(side_length=20.0)

        plane = NumberPlane(
            x_range=(-3, 3),
            y_range=(-0.5, 2),
            x_length=15,
            y_length=6.25,
            background_line_style={
                "stroke_color": WHITE,
                "stroke_width": 0,
                "stroke_opacity": 0.6
            },
            axis_config={"color": BLACK}
        )
        sq.move_to([0, 0, 0])
        curve = plane.plot_parametric_curve(
            lambda t: np.array(
                [
                    joukouski(shiftfunction(circle(t))).real,
                    joukouski(shiftfunction(circle(t))).imag,
                    0,
                ]
            ),
            t_range=[0, 2 * np.pi],
            color="#ece6e2",
            stroke_width=0
        )
        un = Intersection(sq, curve, color="#ece6e2",
                          fill_opacity=1, stroke_color=BLACK)
        self.add(plane, curve, un)


class uniformFlow(Scene):
    def construct(self):
        self.camera.background_color = "#FFF"
        def func(pos): return ((80*RIGHT))
        self.add(ArrowVectorField(func, min_color_scheme_value=u0,
                                  max_color_scheme_value=(u0+10)).scale(2))


class vortex(Scene):
    def construct(self):
        self.camera.background_color = "#FFF"

        def func(pos):
            if np.abs(pos[0]) < 0.2 and np.abs(pos[1]) < 0.2:
                return 0*UR
            else:
                return 10*((-1)*pos[1]/(2*np.pi*(pos[0]**2 + pos[1]**2))*RIGHT + (pos[0]/(2*np.pi*(pos[0]**2 + pos[1]**2)))*UP)
        vf = ArrowVectorField(func, opacity=25, min_color_scheme_value=1.2,
                              stroke_width=7, max_color_scheme_value=3,).scale(2.3)
        self.add(vf)


# edit this to add more genreal cases
class FlowAroundCircle(Scene):
    def construct(self):
        self.camera.background_color = "#FFF"

        def func(pos):
            if np.abs(pos[0]) < 1 and np.abs(pos[1]) < 1:
                return 0*UR
            else:
                return ((2*u0*pos[1]**2/(pos[0]**2 + pos[1]**2)**2 + u0 - u0/(pos[0]**2 + pos[1]**2))*RIGHT + (-2*u0*pos[0]*pos[1]/(pos[0]**2 + pos[1]**2)**2)*UP)/80
        vf = ArrowVectorField(func, stroke_width=25, min_color_scheme_value=((u0-20)/80),
                              max_color_scheme_value=(u0+20)/80).scale(2.3)

        plane = NumberPlane(background_line_style={
            "stroke_color": WHITE,
            "stroke_width": 0,
            "stroke_opacity": 0.6
        },
            axis_config={"color": BLACK})

        stream = StreamLines(
            func, stroke_width=3, max_anchors_per_line=20, min_color_scheme_value=((u0-20)/80),
            max_color_scheme_value=(u0+20)/80
        ).scale(2.3).move_to(
            plane.coords_to_point(3, 0))
        circle = Circle(radius=2.3, color="#ece6e2", stroke_color=BLACK, fill_opacity=1).move_to(
            plane.coords_to_point(0, 0))

        self.add(vf, stream, plane, circle,)


class StreamLinesCircle(Scene):
    def construct(self):
        def func(pos):
            div = 200
            if np.abs(pos[0]) < 1 and np.abs(pos[1]) < 1:
                return 0*UR
            else:
                return ((200*pos[1]**2/(pos[0]**2 + pos[1]**2)**2 + 100 - 100/(pos[0]**2 + pos[1]**2))*RIGHT + (-200*pos[0]*pos[1]/(pos[0]**2 + pos[1]**2)**2)*UP)/div
        vf = ArrowVectorField(func)
        plane = NumberPlane()
        self.add(vf, plane, StreamLines(
            func, stroke_width=3, max_anchors_per_line=30))


class StreamLinesCircleCirculation(Scene):
    def construct(self):
        self.camera.background_color = "#FFF"

        def func(pos):
            div = 100
            if np.abs(pos[0]) < 1 and np.abs(pos[1]) < 1:
                return 0*UR
            else:
                return ((gamma*pos[1]/(2*np.pi*(pos[0]**2 + pos[1]**2)) - 2*R**2*u0*pos[0]**2/(pos[0]**2 + pos[1]**2)**2 + R**2*u0/(pos[0]**2 + pos[1]**2) + u0)*RIGHT + (-gamma*pos[0]/(2*np.pi*(pos[0]**2 + pos[1]**2)) - 2*R**2*u0*pos[0]*pos[1]/(pos[0]**2 + pos[1]**2)**2)*UP)/div
        vf = ArrowVectorField(func, min_color_scheme_value=((u0-15)/100),
                              max_color_scheme_value=((u0+15)/100))
        plane = NumberPlane(
            background_line_style={
                "stroke_color": WHITE,
                "stroke_width": 0,
                "stroke_opacity": 0.6
            },
            axis_config={"color": BLACK}
        )
        circle = Circle(radius=1, color="#ece6e2", stroke_color=BLACK, fill_opacity=1).move_to(
            plane.coords_to_point(0, 0))
        items = VGroup(vf, plane, StreamLines(
            func, stroke_width=3, max_anchors_per_line=30, min_color_scheme_value=((u0-15)/100),
            max_color_scheme_value=((u0+15)/100)), circle)
        self.add(items.scale(2))


def potential1(x, y):
    res = gamma*np.log(np.abs(eta0 - x/2 - complex(0, 1)*y/2 + np.sqrt(-4*a**2 + x**2 + 2*complex(0, 1)*x*y - y**2)/2))/(2*np.pi) + R**2*u0*(-y/2 + (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.sin(np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 + (eta0.imag))/((x/2 - (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.cos(
        np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 - (eta0.real))**2 + (y/2 - (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.sin(np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 - (eta0.imag))**2) + u0*(y/2 - (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.sin(np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 - (eta0.imag))
    return res


def potential2(x, y):
    res = gamma*np.log(np.abs(-eta0 + x/2 + complex(0, 1)*y/2 + np.sqrt(-4*a**2 + x**2 + 2*complex(0, 1)*x*y - y**2)/2))/(2*np.pi) + R**2*u0*(-y/2 - (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.sin(np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 + (eta0.imag))/((x/2 + (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.cos(
        np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 - (eta0.real))**2 + (y/2 + (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.sin(np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 - (eta0.imag))**2) + u0*(y/2 + (4*x**2*y**2 + (-4*a**2 + x**2 - y**2)**2)**(1/4)*np.sin(np.arctan2(2*x*y, -4*a**2 + x**2 - y**2)/2)/2 - (eta0.imag))
    return res


def derivate(func, x0, y0, param):
    h = 0.0001
    if param == "x":
        return (func(x0+h, y0) - func(x0, y0))/h
    else:
        return (func(x0, y0+h) - func(x0, y0))/h


class StreamLinesAirfoil(Scene):
    def construct(self):
        self.camera.background_color = "#FFF"
        sq = Square(side_length=20.0)
        sq.move_to([0, 0, 0])

        def func1(pos):
            x = pos[0]
            y = pos[1]
            p = complex(x, y)
            div = 250
            if np.sqrt(((inverseJoukouski2(inverseShiftfunction(p))).real)**2 + ((inverseJoukouski2(inverseShiftfunction(p))).imag)**2) >= 1.2 and np.sqrt(((inverseJoukouski1(inverseShiftfunction(p))).real)**2 + ((inverseJoukouski1(inverseShiftfunction(p))).imag)**2) >= 1.2:
                return 0 * UR
            else:
                if x < -eta0.real:
                    return ((derivate(potential1, x, y, "y"))*RIGHT + (-1*derivate(potential1, x, y, "x"))*UP)/div
                elif x > -eta0.real and x < eta0.real:
                    return (((derivate(potential2, x+0.3, y+0.3, "y")))*RIGHT + (0)*UP)/div
                else:
                    return ((derivate(potential2, x, y, "y"))*RIGHT + (-1*derivate(potential2, x, y, "x"))*UP)/div

        vf = ArrowVectorField(
            func1,
            min_color_scheme_value=(u0/250),
            max_color_scheme_value=(u0+10)/250
        )
        stream = StreamLines(
            func1, stroke_width=3, max_anchors_per_line=10, min_color_scheme_value=(u0/250),
            max_color_scheme_value=(u0+10)/250
        )
        plane = NumberPlane(
            background_line_style={
                "stroke_color": WHITE,
                "stroke_width": 0,
                "stroke_opacity": 0.6
            },
            axis_config={
                "color": BLACK,
                "stroke_width": 2,
            }
        )
        curve = plane.plot_parametric_curve(
            lambda t: np.array(
                [
                    joukouski(shiftfunction(circle(t))).real,
                    joukouski(shiftfunction(circle(t))).imag,
                    0,
                ]
            ),
            t_range=[0, 2 * np.pi],
            color="#000",
            stroke_width=0
        )
        area = Intersection(sq, curve,              color="#ece6e2",
                            fill_opacity=1, stroke_color=BLACK)
        self.add_foreground_mobjects(area)
        items = VGroup(vf, plane, area, stream)
        items.scale(2)
        self.add(items)


with tempconfig({"quality": "medium_quality", "preview": False, "pixel_width": 1920, "pixel_height": 1080}):
    scene = StreamLinesAirfoil()
    scene.render()
