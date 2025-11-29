#!/bin/bash

compensations=(0 1 2 3)
sharpnesses=(0.001 50)

for comp in "${compensations[@]}"; do
  for sharp in "${sharpnesses[@]}"; do
    echo "➡️  Launching: compensation=$comp, sharpness=$sharp"
    bash sim_setup_v2.sh "$comp" "$sharp"
  done
done
