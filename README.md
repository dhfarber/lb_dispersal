# SLIR Late Blight Model
Model_2D is a Python script to model the spread of a wide range of diseases via a compartmental time-step model.

Model_2d_102818 is the latest version, which can accept non-square fields in addition to square fields.

lb_auto_c25 iteratively runs 27 simulations - a complete factorial of 3 levels each for DMFR (reproduction rate per day), initial infection, and lesion growth rate - and saves them as .pickle files; takes about 24 hours on an average PC.
