#!/bin/bash

#archivo QG_barotrop.f tiene alpha=beta=gama=0, voy cambiando los valores de a uno 
cd ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran

sed -i 's/alpha=0/alpha=2/g' QG_barotrop.f       # alpha=2, beta=0, gama=0
gfortran QG_barotrop.f -o QG
rm ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp/*.dat
./QG
cp -R ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim1

sed -i 's/alpha=2/alpha=0/g; s/beta=0/beta=2/g' QG_barotrop.f    # alpha=0, beta=2, gama=0
gfortran QG_barotrop.f -o QG
rm ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp/*.dat
./QG
cp -R ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim2

sed -i 's/beta=2/beta=0/g; s/gama=0/gama=2/g' QG_barotrop.f      # alpha=0, beta=0, gama=2
gfortran QG_barotrop.f -o QG
rm ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp/*.dat
./QG
cp -R ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim3

sed -i 's/alpha=0/alpha=1/g; s/beta=0/beta=1/g; s/gama=2/gama=0/g' QG_barotrop.f     # alpha=1, beta=1, gama=0
gfortran QG_barotrop.f -o QG
rm ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp/*.dat
./QG
cp -R ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim4

sed -i 's/alpha=1/alpha=0/g; s/gama=0/gama=1/g' QG_barotrop.f    # alpha=0, beta=1, gama=1
gfortran QG_barotrop.f -o QG
rm ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp/*.dat
./QG
cp -R ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim5

sed -i 's/alpha=0/alpha=1/g; s/beta=1/beta=0/g' QG_barotrop.f    # alpha=1, beta=0, gama=1
gfortran QG_barotrop.f -o QG
rm ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp/*.dat
./QG
cp -R ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim6

sed -i 's/alpha=1/alpha=0.66/g; s/beta=0/beta=0.66/g; s/gama=1/gama=0.66/g' QG_barotrop.f    # alpha=0.66, beta=0.66, gama=0.66
gfortran QG_barotrop.f -o QG
rm ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp/*.dat
./QG
cp -R ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran/out_tmp ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/simulaciones/sim7

source activate daniuenv
cd ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/rutinas
python ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/rutinas/plots.py

cd ~/Documentos/materias/metodos_numericos/monografia_final/modelo_QG/modelo_fortran
sed -i 's/alpha=0.66/alpha=0/g; s/beta=0.66/beta=0/g; s/gama=0.66/gama=0/g' QG_barotrop.f    # alpha=0, beta=0, gama=0

