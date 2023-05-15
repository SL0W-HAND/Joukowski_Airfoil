from manim import *
import numpy as np


def circle(theta):
    return np.exp(complex(0, PI*theta))


def shiftfunction(function):
    return complex(0.2, 0.15) - function


def joukouski(z):
    return z + ((1.2)**2)/z


class JoukouskiTransform(Scene):
    def construct(self):
        # create two Cartesian planes
        plane1 = NumberPlane(x_range=(-3, 3), y_range=(-3, 3))
        plane2 = NumberPlane(x_range=(-3, 3), y_range=(-3, 3))
        plane3 = NumberPlane(x_range=(-3, 3), y_range=(-3, 3))

        # position the axes
        plane1.to_corner(UL)
        plane3.to_corner(UR)
        plane2.to_edge(DOWN)
        axes = VGroup(plane1, plane2, plane3)

        self.add(axes)
        # create the arrow
        arrow = Arrow(start=plane1.get_right(),
                      end=plane2.get_left(), buff=SMALL_BUFF)
        label1 = Tex(r"$Z=z_0 - z$", font_size=20).next_to(arrow, DOWN)
        # add the arrow to the scene

        arrow2 = Arrow(start=plane1.get_right(),
                       end=plane3.get_left(), buff=SMALL_BUFF)
        label2 = Tex(
            r"$w = z_0 - z-\frac{a}{z_0 - z}$", font_size=20).next_to(arrow2, UP)

        arrow3 = Arrow(start=plane2.get_right(),
                       end=plane3.get_left(), buff=SMALL_BUFF)
        label3 = Tex(r"$w = Z-\frac{a}{Z}$",
                     font_size=20).next_to(arrow3, DOWN)
        self.add(arrow, label1, arrow2, label2, arrow3, label3)

        curve = plane1.plot_parametric_curve(
            lambda t: np.array(
                [
                    circle(t).real,
                    circle(t).imag,
                    0,
                ]
            ),
            t_range=[0, 2 * PI],
            color="#0FF1CE",
        )

        curve2 = plane2.plot_parametric_curve(
            lambda t: np.array(
                [
                    shiftfunction(circle(t)).real,
                    shiftfunction(circle(t)).imag,
                    0,
                ]
            ),
            t_range=[0, 2 * PI],
            color="#0FF1CE",
        )

        curve3 = plane3.plot_parametric_curve(
            lambda t: np.array(
                [
                    joukouski(shiftfunction(circle(t))).real,
                    joukouski(shiftfunction(circle(t))).imag,
                    0,
                ]
            ),
            t_range=[0, 2 * PI],
            color="#0FF1CE",
        )

        self.add(curve, curve2, curve3)


class Airfoil(Scene):
    def construct(self):

        plane = NumberPlane(x_range=(-3, 3), y_range=(-0.5, 2))
        curve = plane.plot_parametric_curve(
            lambda t: np.array(
                [
                    joukouski(shiftfunction(circle(t))).real,
                    joukouski(shiftfunction(circle(t))).imag,
                    0,
                ]
            ),
            t_range=[0, 2 * PI],
            color="#0FF1CE",
        )
        self.add(plane, curve)


class uniformFlow(Scene):
    def construct(self):
        def func(pos): return ((80*RIGHT))

        self.add(ArrowVectorField(func))


class vortex(Scene):
    def construct(self):
        def func(pos):
            if np.abs(pos[0]) < 0.5 and np.abs(pos[1]) < 0.5:
                return 0*UR
            else:
                return 10*((-1)*pos[1]/(2*np.pi*(pos[0]**2 + pos[1]**2))*RIGHT + (pos[0]/(2*np.pi*(pos[0]**2 + pos[1]**2)))*UP)
        vf = ArrowVectorField(func)
        self.add(vf)


class FlowAroundCircle(Scene):
    def construct(self):
        def func(pos):
            if np.abs(pos[0]) < 1 and np.abs(pos[1]) < 1:
                return 0*UR
            else:
                return (200*pos[1]**2/(pos[0]**2 + pos[1]**2)**2 + 100 - 100/(pos[0]**2 + pos[1]**2))*RIGHT + (-200*pos[0]*pos[1]/(pos[0]**2 + pos[1]**2)**2)*UP
        vf = ArrowVectorField(func)
        plane = NumberPlane()
        self.add(vf, plane)


with tempconfig({"quality": "medium_quality", "preview": False, "pixel_width": 1920, "pixel_height": 1080}):
    scene = FlowAroundCircle()
    scene.render()
