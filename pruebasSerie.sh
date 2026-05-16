#/bin/bash
echo "> T Resultados Serie" >resultadosSerie.log
for a in {1..3}; do
   ssh nodo-0-52 /home/scp06/proyecto/serie/bin.sh >> resultadosSerie.log
 done
