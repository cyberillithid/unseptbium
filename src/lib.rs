use pyo3::prelude::*;

pub fn test_reference(a: usize, b: usize) -> Result<usize, ()> {
    Ok(a+b)
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    let res = test_reference(a, b);
    if let Ok(ret) = res {
        Ok(ret.to_string())
    } else {
        Ok("ERROR OCCURED".to_string())
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn ss13(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    Ok(())
}
