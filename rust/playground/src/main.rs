use std::f64::consts::PI;
use num::complex::Complex;

fn main() {
    println!("Hello, world!");
    
    let mut vector: Vec<i32> = vec![1, 2, 3, 4];
    vector.push(5);
    let slice: &[i32] = &vector;
    println!("{:?} {:?}", vector, slice);
    
    // A tuple is a fixed-size set of values of possibly different types
    let x: (i32, &str, f64) = (1, "hello", 3.4);

    // Destructuring `let`
    let (a, b, c) = x;
    println!("{} {} {}", a, b, c); // 1 hello 3.4
    
    // Indexing
    println!("{}", x.1); // hello


    // Struct
    struct Point {
        x: i32,
        y: i32,
    }

    let origin: Point = Point { x: 10, y: -20 };
    
    println!("{} {}", origin.x, origin.y);
    
    let x = Complex::new(0.0, 2.0*PI);

    println!("e^(2i * pi) = {}", x.exp()); // =~1
    
    let vc: Vec<Complex<f64>> = vec![Complex::new(0.1, 0.2), Complex::new(0.2, 0.3)];
    
    println!("{:?}", vc);
    
    let zero_vec: Vec<Complex<i32>> = vec![Complex::new(0, 0); 5];
    
    println!("{:?}", zero_vec);
    
}
