use std::collections::HashMap;

use ss13::test_reference;
use dmm_tools::dmm::Map;
fn main() {
    let result = test_reference(1400, 88).unwrap();
    println!("Hallo, {}", result);
    let map = Map::from_file(
        std::path::Path::new("E:\\Code\\opensource\\games\\ScavStation\\maps\\ministation\\ministation-0.dmm")).unwrap();
    println!("X,Y,Z: {:?}; len: {} ", map.dim_xyz(), map.key_length());
    // Temporary shit: (118,45)..(133,76)
    let mut ret = HashMap::new();
    for i in 116..135  {
        for j in 41..78 {
                let ix = [0, 255-j, i-1];
                let k = map.grid[ix];
                // println!("Key for SMCore[{:?}]: {:?}", ix,k);
                // println!("\t {:?}", map.dictionary[&k]);
                let prefvecs = map.dictionary[&k].clone();
                let valvec:Vec<(String, HashMap<String, String>)> = prefvecs.into_iter().map(|x| {
                    let hmaps: HashMap<String, String> = x.vars.into_iter().map(|(k, v)| (k, format!("{:?}", v))).collect();
                    (x.path, hmaps)
                }).collect();
                ret.insert(i*256|j, valvec);
            // }
        }
    } // cycle
    let s = serde_json::to_string(&ret).unwrap();
    let _ = std::fs::write("./map.json", s);
}