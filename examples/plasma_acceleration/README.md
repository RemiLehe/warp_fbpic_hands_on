# Laser-plasma examples

## 1. Plasma-wakefield acceleration with FBPIC (in the lab frame)

- Open the file `fbpic_script.py` in
`examples/plasma_acceleration/1_pwfa_labframe_fbpic`
- The different functions/classes used are explained in [the online documentation](https://fbpic.github.io/api_reference/api_reference.html)
- Run the file with
```
python fbpic_script.py
```

**Post-processing:**
- If you run on your local computer, type:
```
jupyter notebook openPMD-visualization.ipynb
```
- If you run on Binder, use the file navigator (not the Terminal)
in order to go to the folder `examples/plasma_acceleration/1_pwfa_labframe_fbpic`
and click on `openPMD-visualization.ipynb`.


**For other, more complete FBPIC examples (e.g. including ionization), see [this page](https://fbpic.github.io/how_to_run.html)**


**For more info on how to use `openPMD-viewer`, see [this page](https://github.com/openPMD/openPMD-viewer/tree/master/tutorials)**

## 2. Laser-wakefield acceleration with Warp (in the lab frame)

- Open the file `warp_script.py` in `examples/plasma_acceleration/2_lwfa_labframe_warp`

**For other, more complete Warp examples, see [this page](https://bitbucket.org/berkeleylab/warp/src/master/scripts/examples/plasma_acceleration/)**

- Run the file with 
```
python warp_script.py
```

- You can change the dimension to "3d" or "circ" in the file, and run the script again

- For post-processing, use the same steps as for FBPIC


## 3. Laser-wakefield acceleration with FBPIC (in the boosted frame)

For a summary of the boosted-frame technique (from a user's perspective), see
[this page](https://fbpic.github.io/advanced/boosted_frame.html)

- Open the file `fbpic_script.py` in
`examples/plasma_acceleration/3_lwfa_boostedframe_fbpic`
- Run the file with 
```
python fbpic_script.py
```
- Post-processing: as usual