use std::cell::RefCell;
use std::rc::Rc;
use std::cell::RefMut;
use std::collections::VecDeque;

use std::f64::consts::PI;

extern crate term_size;

fn scale_x(x: f64, width: f64, xCam: f64, scale: f64) -> f64 { // 0..x..width -> Mandelbrot X scale -2.5..x..1
  ((x / (width / 3.5)) - 1.75) * scale + xCam
}

fn scale_y(y: f64, height: f64, yCam: f64, scale: f64) -> f64 { // 0..y..height -> Mandelbrot Y scale -1..y..1
  (((y / (height / 2.0)) - 1.0)) * scale + yCam
}

//static MAX_ITERATION: u32 = 160;
//static MAX_ITERATION_SCALE: f64 = 255.0 / (MAX_ITERATION as f64);

/*fn calc_pixel(maxIteration: u32, pX: u32, pY: u32, width: u32, height: u32, xCam: f64, yCam: f64, scale: f64) -> u32 {
  let x0: f64 = scale_x(pX as f64, width as f64, xCam, scale);
  let y0: f64 = scale_y(pY as f64, height as f64, yCam, scale);

  let mut x: f64 = 0.0;
  let mut y: f64 = 0.0;

  let mut iteration: u32 = 0;

  let mut xSq = x * x;
  let mut ySq = y * y;

  while xSq + ySq <= 4.0 && iteration < maxIteration {
    /*let xTemp: f64 = x * x - y * y + x0;
    y = 2.0 * x * y + y0;
    x = xTemp;*/

    let xTemp: f64 = xSq - ySq + x0;
    y = x * y;
    y += y;
    y += y0;

    x = xTemp;

    xSq = x * x;
    ySq = y * y;

    iteration += 1;
  }

  iteration
}*/

struct TracingGlobal {
  data: Vec<f64>,
  done: Vec<u8>,
  queue: VecDeque<u32>,

  width: u32,
  actualWidth: u32,
  height: u32,

  maxIteration: f64,
  periodicityCheckingDelta: f64,

  xCam: f64,
  yCam: f64,
  scale: f64,
  xOffset: u32,

  juliaVectorX: f64,
  juliaVectorY: f64,

  useSmoothColoring: bool
}

static LOG10_2: f64 = 0.301029995663981198017467022509663365781307220458984375; // does not support log10(n) at compile time

static PERIODIC_CHECK_DELTA_DIVIDE: f64 = 10000.0;

fn calc_pixel(tracingGlobal: &TracingGlobal, sX: f64, sY: f64) -> f64 { //pX: u32, pY: u32) -> f64 {
  let mut mandelbrotXOffset: f64 = 0.0;

  let mut x0 = sX;
  let y0 = sY;

  //let mut x0: f64 = scale_x(pX as f64, tracingGlobal.actualWidth as f64, tracingGlobal.xCam, tracingGlobal.scale);
  //let y0: f64 = scale_y(pY as f64, tracingGlobal.height as f64, tracingGlobal.yCam, tracingGlobal.scale);

  let mut juliaX = tracingGlobal.juliaVectorX;
  let mut juliaY = tracingGlobal.juliaVectorY;
  
  if juliaX == 0.0 && juliaY == 0.0 {
    x0 -= 0.5;

    juliaX = x0;
    juliaY = y0;
  }

  let mut x: f64 = x0;
  let mut y: f64 = y0;

  let mut iteration: f64 = 0.0;

  let mut xSq = x * x;
  let mut ySq = y * y;

  let mut r = 0;
  let mut g = 0;
  let mut b = 0;

  let mut x_old: f64 = 0.0;
  let mut y_old: f64 = 0.0;

  let mut period: u32 = 0;

  while xSq + ySq <= 4.0 && iteration < tracingGlobal.maxIteration {
    let xTemp: f64 = xSq - ySq + juliaX;
    y = x * y;
    y += y + juliaY;

    x = xTemp;

    xSq = x * x;
    ySq = y * y;

    iteration += 1.0;

    if (if x > x_old { x - x_old } else { x_old - x }) < tracingGlobal.periodicityCheckingDelta && (if y > y_old { y - y_old } else { y_old - y }) < tracingGlobal.periodicityCheckingDelta {
      iteration = tracingGlobal.maxIteration;
      break;
    }

    period += 1;
    if period > 20 {
      period = 0;

      x_old = x;
      y_old = y;
    }
  }

  if tracingGlobal.useSmoothColoring && iteration < tracingGlobal.maxIteration {
    let log_zn = (xSq + ySq).log10() / 2.0;
    let nu = (log_zn / LOG10_2).log10() / LOG10_2;
    iteration = iteration + 1.0 - nu;
  }

  iteration
}

static DONE_QUEUED: u8 = 2;
static DONE_LOADED: u8 = 1;

fn addQueue(mut tracingGlobal: &mut TracingGlobal, p: usize) {
  if tracingGlobal.done[p] & DONE_QUEUED == DONE_QUEUED {
    return;
  }

  tracingGlobal.done[p] |= DONE_QUEUED;
  tracingGlobal.queue.push_back(p as u32);
}

fn get_xy_from_p(tracingGlobal: &TracingGlobal, p: u32) -> (u32, u32) {
  (p % tracingGlobal.width, p / tracingGlobal.width)
}

fn mirror_value(val: f64, origin: f64) -> f64 {
  if val > origin { origin - (val - origin) } else { origin + (origin - val) }
}

fn load(tracingGlobal: &mut TracingGlobal, p: usize) -> f64 {
  if tracingGlobal.done[p] & DONE_LOADED == DONE_LOADED {
    return tracingGlobal.data[p];
  }

  let (x, y) = get_xy_from_p(tracingGlobal, p as u32);

  let sX: f64 = scale_x((x + tracingGlobal.xOffset) as f64, tracingGlobal.actualWidth as f64, tracingGlobal.xCam, tracingGlobal.scale);
  let sY: f64 = scale_y(y as f64, tracingGlobal.height as f64, tracingGlobal.yCam, tracingGlobal.scale);

  let result = calc_pixel(tracingGlobal, sX, sY); //calc_pixel(tracingGlobal, x + tracingGlobal.xOffset, y);

  tracingGlobal.done[p] |= DONE_LOADED;

  tracingGlobal.data[p] = result;

  /*let mirrorP = ((mirror_value(y, tracingGlobal.height / 2) * tracingGlobal.width) + x) as usize;*/

  let mirrorSY = mirror_value(sY, 0.0);
  let mirrorPY: u32 = ((0.5 * ((tracingGlobal.height as f64) * (tracingGlobal.scale - tracingGlobal.yCam + mirrorSY))) / tracingGlobal.scale).floor() as u32;

  if mirrorPY < tracingGlobal.height {
    let mirrorP = ((mirrorPY * tracingGlobal.width) + x) as usize;

    tracingGlobal.done[mirrorP] |= DONE_LOADED;

    tracingGlobal.data[mirrorP] = result;
  }

  /*let mirrorSX = mirror_value(sX, 0.125);
  let mirrorPX: u32 = (((tracingGlobal.actualWidth as f64) * ((0.5 * tracingGlobal.scale) - (0.285714 * tracingGlobal.xCam) + (0.285714 * mirrorSX))) / tracingGlobal.scale).floor() as u32 - tracingGlobal.xOffset;

  if mirrorPX < tracingGlobal.width {
    let mirrorP = ((y * tracingGlobal.width) + mirrorPX) as usize;

    tracingGlobal.done[mirrorP] |= DONE_LOADED;

    tracingGlobal.data[mirrorP] = result;
  }*/

  result
}

fn scan(tracingGlobal: &mut TracingGlobal, p: u32) -> f64 {
  let (x, y) = get_xy_from_p(tracingGlobal, p as u32);

  let center = load(tracingGlobal, p as usize);

  let ll = x >= 1;
  let rr = x < tracingGlobal.width - 1;
  let uu = y >= 1;
  let dd = y < tracingGlobal.height - 1;

  let lp = (p - 1) as usize;
  let rp = (p + 1) as usize;
  let up = (p - tracingGlobal.width) as usize;
  let dp = (p + tracingGlobal.width) as usize;

  let l = ll && load(tracingGlobal, lp) != center; // tracingGlobal.maxIteration; //center;
  let r = rr && load(tracingGlobal, rp) != center; // tracingGlobal.maxIteration; //center;
  let u = uu && load(tracingGlobal, up) != center; // tracingGlobal.maxIteration; //center;
  let d = dd && load(tracingGlobal, dp) != center; // tracingGlobal.maxIteration; //center;

  if l { addQueue(tracingGlobal, lp); }
  if r { addQueue(tracingGlobal, rp); }
  if u { addQueue(tracingGlobal, up); }
  if d { addQueue(tracingGlobal, dp); }

  if (uu && ll) && (l || u) { addQueue(tracingGlobal, up - 1); }
  if (uu && rr) && (r || u) { addQueue(tracingGlobal, up + 1); }
  if (dd && ll) && (l || d) { addQueue(tracingGlobal, dp - 1); }
  if (dd && rr) && (r || d) { addQueue(tracingGlobal, dp + 1); }

  center
}

fn gen_color_from_iteration(iteration: f64, maxIterationColorScale: f64) -> f64 {
  //float regulator = (iteration as f64) - log(log(dot(p, p)) / log(2.0)) / log(2.0);
  //color = vec3(0.95 + .012 * regulator, 1.0, .1 + .4 * (1.0 + sin(.3 * regulator)));

  iteration * maxIterationColorScale
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
  return (1.0 - t) * a + t * b;
}

fn smooth_color(useSmoothColoring: bool, iteration: f64, maxIterationColorScale: f64) -> u8 {
  let color1 = gen_color_from_iteration(iteration, maxIterationColorScale);

  if useSmoothColoring == false {
    return color1.round() as u8;
  }

  let color2 = gen_color_from_iteration(iteration + 1.0, maxIterationColorScale);

  //color1.round() as u8
  lerp(color1, color2, iteration % 1.0).round() as u8
  //lerp(color1, color2, ((iteration % 5) as f64) / 4.0).round() as u8
}

fn histo_smooth_color(useSmoothColoring: bool, iteration: f64, color1: u8, color2: u8) -> u8 {
  lerp(color1 as f64, color2 as f64, iteration % 1.0).round() as u8
}

pub fn gen_histo_data(maxIteration: f64, juliaVectorX: f64, juliaVectorY: f64, width: u32, height: u32) -> Vec<u8> {
  let areaU: u32 = (width * height);
  let area: usize = areaU as usize;

  let scale = 0.8;

  let mut tracingGlobal = TracingGlobal {
    data: vec![0.0; area],
    done: vec![0; area],

    queue: VecDeque::new(),

    width: width,
    actualWidth: width,
    height: height,

    maxIteration: maxIteration - 1.0,
    periodicityCheckingDelta: scale / PERIODIC_CHECK_DELTA_DIVIDE,

    xCam: 0.0,
    yCam: 0.0,
    scale: scale,
    xOffset: 0,
    
    juliaVectorX: juliaVectorX,
    juliaVectorY: juliaVectorY,

    useSmoothColoring: false
  };

  let refcell: Rc<RefCell<TracingGlobal>> = Rc::new(RefCell::new(tracingGlobal));

  let mut borrowMut = refcell.borrow_mut();

  for y in 0..height {
    let left = y * width;
    addQueue(&mut *borrowMut, left as usize);
    addQueue(&mut *borrowMut, (left + (width - 1)) as usize);
  }

  for x in 1..(width - 1) {
    addQueue(&mut *borrowMut, x as usize);
    addQueue(&mut *borrowMut, ((height - 1) * width + x) as usize)
  }

  let maxIterationU = maxIteration as usize;

  let mut numIterationsPerPixel: Vec<f64> = vec![0.0; maxIterationU];
  let mut iterationsTotal: f64 = 0.0;

  let mut histoColors: Vec<u8> = vec![0; maxIterationU];

  while borrowMut.queue.len() > 0 {
    let p = borrowMut.queue.pop_front().unwrap();

    let it = scan(&mut *borrowMut, p);

    if it != maxIteration - 1.0 {
      numIterationsPerPixel[it as usize] += 1.0;
      iterationsTotal += 1.0;
    }
  }

  for it in 0..maxIterationU {
    let mut out: f64 = 0.0;

    for i in 0..(it + 1) {
      out += numIterationsPerPixel[it] / iterationsTotal;
    }

    histoColors[it] = (out * 255.0) as u8;
  }

  /*for p in 0..area {
    let mut iteration = borrowMut.data[p];

    if p != area - 1 && (borrowMut.done[p] & DONE_LOADED) == 1 {
      if (borrowMut.done[p + 1] & DONE_LOADED) == 0 {
        borrowMut.data[p + 1] = iteration;
        borrowMut.done[p + 1] |= DONE_LOADED;
      }
    }

    if iteration != maxIteration - 1.0 {
      numIterationsPerPixel[iteration as usize] += 1.0;
      iterationsTotal += 1.0;
    }
  }

  for p in 0..area {
    let mut iteration = borrowMut.data[p];
    let iterationU = iteration as usize;

    if histoColors[iterationU] == 0 {
      let mut out: f64 = 0.0;

      for i in 0..(iterationU + 1) {
        out += numIterationsPerPixel[i] / iterationsTotal;
      }

      histoColors[iterationU] = (out * 250.0) as u8;

      //histoColors[iterationU] = (((numIterationsPerPixel[iterationU] / iterationsTotal) * (iteration + 1.0)) * 250.0) as u8;
    }
  }*/

  histoColors//.into_iter().map(JsValue::from).collect()
}

pub fn gen_data(histoColors: &Vec<u8>, linesBetweenColumns: bool, usehistogramColoring: bool, useSmoothColoring: bool, maxIteration: f64, maxIterationColorScale: f64, juliaVectorX: f64, juliaVectorY: f64, renderWidth: u32, actualWidth: u32, xOffset: u32, height: u32, xCam: f64, yCam: f64, scale: f64) -> Vec<u8> {
  //set_panic_hook();

  //log(&smooth_color(10, maxIterationColorScale).to_string());

  let areaU: u32 = (renderWidth * height);
  let area: usize = areaU as usize;

  let mut pixels: Vec<u8> = vec![0; area * 3]; //Uint8ClampedArray = Uint8ClampedArray::new_with_length(areaU * 4);

  let mut tracingGlobal = TracingGlobal {
    data: vec![0.0; area],
    done: vec![0; area],

    queue: VecDeque::new(),

    width: renderWidth,
    actualWidth: actualWidth,
    height: height,

    maxIteration: maxIteration - 1.0,
    periodicityCheckingDelta: scale / PERIODIC_CHECK_DELTA_DIVIDE,

    xCam: xCam,
    yCam: yCam,
    scale: scale,
    xOffset: xOffset,
    
    juliaVectorX: juliaVectorX,
    juliaVectorY: juliaVectorY,

    useSmoothColoring: useSmoothColoring
  };

  let mut i = 0;

  /*if !useBorderTracing {
    for y in 0..height {
      for x in 0..renderWidth {
        let iteration = calc_pixel(&tracingGlobal, x + xOffset, y); //calc_pixel(maxIteration, x + xOffset, y, actualWidth, height, xCam, yCam, scale);
  
        let mut r = 0;
        let mut g = 0;
        let mut b = 0;
  
        if iteration != maxIteration {
          let color = gen_color_from_iteration(iteration, maxIterationColorScale);
          r = color;
          g = color;
          b = 255;
        }
  
        pixels.set_index(i, r);
        pixels.set_index(i + 1, g);
        pixels.set_index(i + 2, b);
        pixels.set_index(i + 3, if linesBetweenColumns && x == renderWidth - 1 { 0 } else { 255 });
  
        i += 4;
      }
    }

    return pixels;
  }*/
  

  let refcell: Rc<RefCell<TracingGlobal>> = Rc::new(RefCell::new(tracingGlobal));

  let mut borrowMut = refcell.borrow_mut();

  for y in 0..height {
    let left = y * renderWidth;
    addQueue(&mut *borrowMut, left as usize);
    addQueue(&mut *borrowMut, (left + (renderWidth - 1)) as usize);
  }

  for x in 1..(renderWidth - 1) {
    addQueue(&mut *borrowMut, x as usize);
    addQueue(&mut *borrowMut, ((height - 1) * renderWidth + x) as usize)
  }

  let maxIterationU = maxIteration as usize;

  //let mut numIterationsPerPixel: Vec<f64> = vec![0.0; maxIterationU];
  //let mut iterationsTotal: f64 = 0.0;

  while borrowMut.queue.len() > 0 {
    let p = borrowMut.queue.pop_front().unwrap();

    let it = scan(&mut *borrowMut, p);

    /*if usehistogramColoring && it != maxIteration - 1.0 {
      numIterationsPerPixel[it as usize] += 1.0;
      iterationsTotal += 1.0;
    }*/
  }


  //let mut histoColors: Vec<u8> = vec![0; maxIterationU];

  for p in 0..area {
    let mut iteration = borrowMut.data[p];

    /*if usehistogramColoring && iteration < 25 {
      //iteration = if maxIteration < 20 { maxIteration - 1 } else { 10 };
      histoColors[iterationU] = 200;
    }*/

    if p != area - 1 && (borrowMut.done[p] & DONE_LOADED) == 1 {
      if (borrowMut.done[p + 1] & DONE_LOADED) == 0 {
        borrowMut.data[p + 1] = iteration;
        borrowMut.done[p + 1] |= DONE_LOADED;
      }
    }

    let mut r = 0;
    let mut g = 0;
    let mut b = 0;
    let mut a = if linesBetweenColumns && ((p as u32) % renderWidth) == renderWidth - 1 { 0 } else { 255 };

    if usehistogramColoring {
      /*if histoColors[iterationU] == 0 {
        if iteration != maxIteration - 1.0 {
          histoColors[iterationU] = ((((numIterationsPerPixel[iterationU] / iterationsTotal) * (iteration + 1.0)) * 255.0) as u8);
        } else {
          histoColors[iterationU] = 255;
        }
      }*/

      if iteration < 10.0 {
        iteration = 10.0;
      }

      let iterationU = iteration.floor() as usize;

      if iteration != maxIteration - 1.0 {
        let color = histo_smooth_color(useSmoothColoring, iteration, histoColors[iterationU], histoColors[iterationU + 1]);

        r = 255 - color;
        g = 255 - color;
        b = 255;
      }
    } else {
      if iteration != maxIteration - 1.0 {
        let color = smooth_color(useSmoothColoring, iteration, maxIterationColorScale); //gen_color_from_iteration(iteration, maxIterationColorScale);
        r = color;
        g = color;
        b = 255;
      }
    }

    pixels[i] = r;
    pixels[i + 1] = g;
    pixels[i + 2] = b;
    //pixels[i + 3]

    i += 3;
  }

  pixels
}

pub fn main() {
  let maxIteration = 3000.0;
  let maxIterationColorScale = 255.0 / maxIteration;

  let juliaVectorX = 0.0;
  let juliaVectorY = 0.0;

  let (_width, _height) = term_size::dimensions().unwrap();

  let width = _width as u32;
  let height = _height as u32;

  let mut xCam = 0.0;
  let yCam = 0.0;

  let mut scale = 1.2;

  let histoData = gen_histo_data(maxIteration, juliaVectorX, juliaVectorY, width, height * 2);

  let refcell: Rc<RefCell<Vec<u8>>> = Rc::new(RefCell::new(histoData));

  let mut frame = 0;

  while (true) {
    let mut borrowMut = refcell.borrow();

    xCam -= 0.005;
    scale *= 0.98;

    let data = gen_data(& *borrowMut, false, true, false, maxIteration, maxIterationColorScale, juliaVectorX, juliaVectorY, width, width, 0, height * 2, xCam, yCam, scale);

    let mut lastYValues: Vec<u8> = vec![254; _width * 3];

    let printChar = &"▄";

    let mut y = 0;

    let area = _width * (_height * 2);

    for j in 0..area {
      let i = j * 3;
      let r = data[i];
      let g = data[i + 1];
      let b = data[i + 2];

      if j % _width == 0 {
        y += 1;

        if y > 1 {
          y = 0;
        }
      }

      let x = j % _width * 3;

      if y == 0 {
        lastYValues[x] = r;
        lastYValues[x + 1] = g;
        lastYValues[x + 2] = b;

        continue;
      }

      let mut lastR = lastYValues[x];
      let mut lastG = lastYValues[x + 1];
      let mut lastB = lastYValues[x + 2];

      if lastR == 254 {
        lastR = r;
        lastG = g;
        lastB = b;
      }

      if j > area - 3 {
        continue;
      }

      print!("{}", background_rgb(lastR, lastG, lastB, &foreground_rgb(r, g, b, printChar)));
    }

    print!("{esc}[1;1H", esc = 27 as char);

    frame += 1;
  }
}

fn foreground_rgb(r: u8, g: u8, b: u8, text: &str) -> String {
  format!("\x1b[38;2;{};{};{}m{}\x1b[0m", r, g, b, text)
}

fn background_rgb(r: u8, g: u8, b: u8, text: &str) -> String {
  format!("\x1b[48;2;{};{};{}m{}\x1b[0m", r, g, b, text)
}