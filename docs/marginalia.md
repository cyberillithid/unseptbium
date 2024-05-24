# Math model of supermatter core engine

useless contour:
- full 16x7 + 1
- perim 17+14+14
= total 158 cells

cold contour:
- in 10
- full 15x5 = 75
- full 11x2 = 22 
- out 1 + 8 + 5 = 14
= total 121 cells

Contour 1 / cold:
1. right N2 canister -> pump -> green 16 cells
2. join with cold-contour-'out' + 1
3. high_power pump + 2 -> cold turbine -> black + 16
4. cold-contour-'in' -> radiators cold -> +5 (2)

Contour 2 / hot:
1. left N2 canister -> pump -> blue (10 cells)
2. join with hotloop -> pipe 8 cells
3.1 6 pipes to scrubbers -> waste
3.2 core injector -> core
4. core -> core vent -> 4 pipes -> hot turbine -> 4 pipes -> (2)



-----

sidenotes

`[X]` -- molar concentration, in moles/cm3 -- well describes chem reactions
that's exactly what EPR shows 


32x48 -- supermatter sprite
so we deem it to be a bit less than full 2.5 x 1x1 m ?


ref.stat -- expansion... into hot loop, I think. or something. yeah, bad choice. gotta re-`debug` it for cold loop later.


## Things to change (probably)

/datum/codex_entry/supermatter

SM core monitoring template to use `Ex*` instead of `PH`
@ `nano\templates\supermatter_monitor.tmpl` : 63


metodika-teplovogo-raschyota-teploobmennogo-apparata-dlya-ohlazhdeniya-vozduha-pered-kompressorom-dvigatelya-vozdushno-kosmicheskogo-samolyota.pdf
