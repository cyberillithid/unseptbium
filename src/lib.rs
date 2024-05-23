use std::collections::HashMap;

use dreammaker::{constants::Constant, DMError};
use pyo3::{exceptions::PyTypeError, prelude::*}; //  types::{PyList, PyNone, PyTuple}

use dmm_tools::dmm::Map; // , Prefab

pub fn test_reference(a: usize, b: usize) -> Result<usize, ()> {
    Ok(a+b)
}

#[pyclass]
/// Thin wrapper over SpaceManiac's DMM tools Map Rust class
struct Ss13Map {
    inner: Result<Map, DMError>,
    path: String
}

#[pyclass]
#[derive(PartialEq)]
/// A minimal sub-implementation, because I don't need finer details for now.
enum ByondConst {
    Null { },
    String { s : String },
    Resource { rsrc : String },
    Float { f : f32 },
    // not actually implemented
    List { repr: String },
    New { arg1 : String, arg2: String },
    Call { arg1 : String, arg2: String },
    Prefab { arg1 : String, arg2: String },
}

impl From<Constant> for ByondConst {
    fn from(value: Constant) -> Self {
        match value {
            Constant::Null(_) => ByondConst::Null {  },
            Constant::Float(f) => ByondConst::Float { f },
            Constant::String(s) => ByondConst::String { s: s.into() },
            Constant::Resource(s) => ByondConst::Resource{ rsrc: s.into() },
            // not implemented yet
            Constant::New { type_, args } => 
                ByondConst::New { 
                    arg1: format!("{:?}", type_), 
                    arg2: format!("{:?}", args) }
            ,
            Constant::List(l) => ByondConst::List { 
                repr: format!("{:?}", l) },
            Constant::Call(f, args) => 
                ByondConst::Call{
                 arg1: f.to_string(), arg2: format!("{:?}", args) },
            Constant::Prefab(p) => ByondConst::Prefab { arg1: format!("{:?}", p.path),
             arg2: format!("{:?}", p.vars) },
        }
    }
}

#[pymethods]
impl ByondConst {
    fn __repr__(&self) -> String {
        match self {
            ByondConst::Null {  } => "ByondConst()".into(),
            ByondConst::String { s } => format!("\"{}\"", s),
            ByondConst::Resource { rsrc } => format!("ByondConst.Resource({})", rsrc),
            ByondConst::Float { f } => format!("{}", f),
            ByondConst::List { repr } => format!("ByondConst.List({})", repr),
            ByondConst::New { arg1, arg2 } => format!("ByondConst.New({}, {})", arg1, arg2),
            ByondConst::Call { arg1, arg2 } => format!("ByondConst.Call({}, {})", arg1, arg2),
            ByondConst::Prefab { arg1, arg2 } => format!("ByondConst.Prefab({}, {})", arg1, arg2),
        }
    }
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
    fn __int__ (&self) -> PyResult<isize> {
        match self {
            ByondConst::Null {  } => Ok(0),
            ByondConst::Float { f } => Ok(f.round() as isize),
            _ => Err(PyErr::new::<PyTypeError, _>("Not a numeric"))
        }
    }
}

#[pymethods]
impl Ss13Map {
    #[new]
    fn new(path: String) -> Self {
        Self {
            inner: Map::from_file(
                std::path::Path::new(path.as_str())),
            path
        }
    }

    fn __repr__(&self) -> String {
        format!("Ss13Map({})", self.path)
    }

    fn __str__(&self) -> String{
        if self.inner.is_ok() {
            format!("map {}", self.path)
        } else {
            format!("map LOAD FAILED {}", self.path)
        }
    }

    fn __getitem__(&self, a: (usize,usize)) -> Vec<(String, HashMap<String, ByondConst>)> {
        if let Ok (pmap) = self.inner.as_ref() {
            let (i,j) = a;
            let grid_ix = [0, 255-j, i-1];
            let k = pmap.grid[grid_ix];
            let prefabs = pmap.dictionary[&k].clone();
            prefabs.into_iter().map(
                |x| (x.path, x.vars.into_iter().map(
                    |(k,v)| (k,v.into())
                ).collect())).collect()
            // Vec::new()
        } else {
            Vec::new()
        }
    }

    /// Finds all indexes which contain given prefab with given path.
    fn find(&self, path: String) -> Vec<(usize, usize)> {
        if let Ok (pmap) = self.inner.as_ref() {
            // F
            let ks: Vec<_> = pmap.dictionary.iter().filter_map(
                |(k,v)| 
                    v.iter().find(|pf| pf.path == path).map(|_| k.clone())
                ).collect();
            pmap.grid.indexed_iter().filter_map(
                |((_,j,i), v)| 
                    ks.iter().find(|x| *x==v).map(|_| (i+1, 255-j))
            ).collect()
        } else {
            Vec::new()
        }
    }
}

/// Formats the sum of two numbers as string.
// #[pyfunction]
// fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
//     let res = test_reference(a, b);
//     if let Ok(ret) = res {
//         Ok(ret.to_string())
//     } else {
//         Ok("ERROR OCCURED".to_string())
//     }
// }

/// A Python module for SS13 .dmm parsing 
/// (mainly for pneumohydrothermal schematics purposes)
/// implemented in Rust.
#[pymodule]
fn ss13_usb (m: &Bound<'_, PyModule>) -> PyResult<()> {
    // m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<Ss13Map>()?;
    m.add_class::<ByondConst>()?;
    Ok(())
}
