package main

import (
	"errors"
	"fmt"
	"github.com/faiface/pixel"
	"github.com/faiface/pixel/imdraw"
	"github.com/faiface/pixel/pixelgl"
	"math"
	"math/cmplx"
)

var width,height = 750.0, 750.0
var xmin,xmax = -2.0,1.0
var ymin,ymax = -1.5,1.5
const s = 1
var imd = imdraw.New(nil) // Not using textured polygons.

func mapToComplex(x,y float64) (complex128) {
	return complex(xmin + (float64(x)/float64(width))*(xmax-xmin), ymax - (float64(y)/float64(height))*(ymax-ymin))
}

func f(z,c complex128) (complex128) {
	return z*z + c
}


// ASSUMES XS IS SORTED
func monoCubicInterpolate(xs []float64, ys []float64) func(float64) float64 {
	length := len(xs)

	// Dealing with length issues
	switch {
	case length != len(ys):
		panic(errors.New("need an equal count of xs and ys"))
	case length == 0:
		return func(x float64) (float64) {
			return 0
		}
	case length == 1:
		return func(x float64) (float64) {
			return ys[0]
		}
	}

	// Get consecutive differences and slopes
	dys := make([]float64, length-1)
	dxs := make([]float64, length-1)
	ms := make([]float64, length-1)
	for i := 0; i < length-1; i++ {
		dx := xs[i+1] - xs[i]
		dy := ys[i+1] - ys[i]
		dxs[i], dys[i], ms[i] = dx, dy, dy/dx
	}

	// Get degree-1 coefficients
	c1s := make([]float64, length)
	c1s[0] = ms[0]
	for i := 1; i < len(dxs); i++ {
		m, mNext := ms[i-1], ms[i]
		if m * mNext <= 0 {
			c1s[i] = 0
		} else {
			dx_, dxNext := dxs[i-1], dxs[i]
			common := dx_ + dxNext
			c1s[i] = 3 * common / ((common + dxNext)/m + (common + dx_)/mNext)
		}
	}
	c1s[length-1] = ms[len(ms)-1]

	// Get degree-2 and degree-3 coefficients
	c2s := make([]float64, length-1)
	c3s := make([]float64, length-1)
	for i := 0; i < len(c1s)-1; i++ {
		c1, m_ := c1s[i], ms[i]
		invDx, common_ := 1/dxs[i], c1 + c1s[i+1] - m_ - m_
		c2s[i] = (m_ - c1 - common_) * invDx
		c3s[i] = common_ * invDx * invDx
	}

	// Return interpolant function
	return func(x float64) float64 {
		// The rightmost point in the dataset should give an exact result
		i := len(xs) - 1
		if x == xs[i] {
			return ys[i]
		}

		// Search for the interval x is in, returning the corresponding y if x is one of the original xs
		low, mid, high := 0, len(c3s)-1, len(c3s)-1
		for low <= high {
			mid = int(math.Floor(0.5 * float64(low + high)))
			xHere := xs[mid]
			if xHere < x {
				low = mid + 1
			} else if xHere > x {
				high = mid - 1
			} else {
				return ys[mid]
			}
		}
		i = int(math.Max(0, float64(high)))

		// Interpolate
		diff := x - xs[i]
		diffSq := diff*diff
		return ys[i] + c1s[i]*diff + c2s[i]*diffSq + c3s[i]*diff*diffSq
	}

}

func getCubicColour(i float64) pixel.RGBA {
	i = math.Mod(i, 1)
	r := monoCubicInterpolate([]float64{0.0, 0.16, 0.42, 0.6425, 0.8575,1.0}, []float64{0, 0.1255, 0.9294, 1.0, 0,0})
	g := monoCubicInterpolate([]float64{0.0, 0.16, 0.42, 0.6425, 0.8575,1.0}, []float64{0.0275, 0.4196, 1.0, 0.6667, 0.0078, 0.0275})
	b := monoCubicInterpolate([]float64{0.0, 0.16, 0.42, 0.6425, 0.8575,1.0}, []float64{0.3922, 0.7961, 1.0, 0, 0, 0.3922})
	return pixel.RGB(r(i), g(i), b(i))
}

func getColour(i float64) (colour pixel.RGBA) {
	var colourAdd = 6.0/100.0
	colVal := math.Mod(i * colourAdd, 6)

	colour.R = math.Min(math.Max(math.Abs(3 - colVal) - 1, 0), 1)
	colour.G = math.Min(math.Max(2 - math.Abs(2 - colVal), 0), 1)
	colour.B = math.Min(math.Max(2 - math.Abs(4 - colVal), 0), 1)
	colour.A = 1
	return
}

func update(win *pixelgl.Window, spacing float64) {

	maxiter :=  math.Floor(math.Sqrt(2 * math.Sqrt(math.Abs(1-math.Sqrt(5 * (1/(xmax-xmin))))))*66.5)

	colours := make([]pixel.RGBA, int(maxiter))
	for i := 0; i < len(colours); i++ {
		colours[i] = getCubicColour(math.Sqrt(float64(i)))
	}

	for py := 0.0; py < height; py+=spacing {
		for px := 0.0; px < width; px+=spacing {
			z := 0+0i
			c := mapToComplex(px, py)

			iter := maxiter
			for i := 0.0; i < maxiter; i++ {
				z = f(z,c)
				// The << operator denotes a bit-shift: a << b = a * 2^b.
				// 2^8 has been chosen as a reasonable bailout radius.
				if cmplx.Abs(z) > (1 << 8) {
					iter = i
					break
				}
			}
			// Used to void floating point issues with points inside the set
			if iter < maxiter {
				// sqrt of inner term removed using log simplification rules.
				log_zn := math.Log(cmplx.Abs(z))
				nu := math.Log(log_zn / math.Log(2)) / math.Log(2)
				iter += 1 - nu
			}

			if iter == maxiter {
				imd.Color = pixel.RGB(0, 0, 0)
			} else {
				imd.Color = getCubicColour(math.Pow(iter/20,0.75))
			}

			imd.Push(pixel.V(px, py))
			imd.Push(pixel.V(px+spacing, py+spacing))
			imd.Rectangle(0)
		}
	}
	imd.Draw(win)
}

func run() {
	cfg := pixelgl.WindowConfig{
		Title: "Mandelbrot Set",
		// Drawing a rectangle within which everything else will be drawn.
		Bounds: pixel.R(0, 0, width, height),
		VSync: true,
	}
	// Creating an actual window
	win, err := pixelgl.NewWindow(cfg)
	if err != nil {
		panic(err)
	}

	update(win, s)

	zoomState := 0
	var coord1, coord2 pixel.Vec
	for !win.Closed() {
		if win.JustPressed(pixelgl.MouseButtonLeft) {
			mouse := win.MousePosition()
			imd.Color = pixel.RGB(1,1,1)
			if zoomState == 0 {
				zoomState += 1
				imd.Push(mouse)
				imd.Push(pixel.V(mouse.X+3,mouse.Y+3))
				imd.Rectangle(0)
				imd.Draw(win)
				coord1 = mouse
			} else if zoomState == 1 {
				zoomState = 0
				imd.Push(coord1)
				distX, distY := math.Abs(mouse.X-coord1.X), math.Abs(mouse.Y-coord1.Y)
				if distX > distY {
					coord2.X = mouse.X
					if mouse.Y > coord1.Y {
						coord2.Y = coord1.Y + distX
					} else {
						coord2.Y = coord1.Y - distX
					}
				} else {
					if mouse.X > coord1.X {
						coord2.X = coord1.X + distY
					} else {
						coord2.X = coord1.X - distY
					}
					coord2.Y = mouse.Y
				}
				imd.Push(coord2)
				imd.Rectangle(3)
				imd.Draw(win)
				win.Update()

				x1 := real(mapToComplex(coord1.X, coord1.Y))
				x2 := real(mapToComplex(coord2.X, coord2.Y))
				y1 := imag(mapToComplex(coord1.X, coord1.Y))
				y2 := imag(mapToComplex(coord2.X, coord2.Y))

				xmin,xmax = math.Min(x1,x2), math.Max(x1,x2)
				ymin,ymax = math.Min(y1,y2), math.Max(y1,y2)
				update(win, s)
			}
		}
		win.Update()
	}
}

func main() {
	var f = func(x int) int {return x*x}
	fmt.Println(f(2))

	pixelgl.Run(run)
}
