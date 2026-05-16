#!/bin/bash
echo "> T Resultados Fase1" > resultadosFase1.log
echo "> T 16 procesos" >> resultadosFase1.log
for a in {1..3}; do
   exec.sh 16 ./paralelo/intentoDeMejoraVersionIvan/bin card >> resultadosFase1.log
done
echo "> T 32 procesos" >> resultadosFase1.log
for a in {1..3}; do
   exec.sh 32 ./paralelo/intentoDeMejoraVersionIvan/bin card >> resultadosFase1.log
done
echo "> T 64 procesos" >> resultadosFase1.log
for a in {1..3}; do
   exec.sh 64 ./paralelo/intentoDeMejoraVersionIvan/bin card >> resultadosFase1.log
done
echo "> T 80" >> resultadosFase1.log
for a in {1..3}; do
   exec.sh 80 ./paralelo/intentoDeMejoraVersionIvan/bin card >> resultadosFase1.log
done
echo "> T 120" >> resultadosFase1.log
for a in {1..3}; do
   exec.sh 120 ./paralelo/intentoDeMejoraVersionIvan/bin card >> resultadosFase1.log
done

echo "> T Resultados Fase2 Jon" > resultadosFase2Jon.log
echo "> T 17 procesos grupos de 4" >> resultadosFase2Jon.log
for a in {1..3}; do
   exec.sh 17 ./paralelo/fase2/bin card 4 >> resultadosFase2Jon.log
done
echo "> T 17 procesos grupos de 8" >> resultadosFase2Jon.log
for a in {1..3}; do
   exec.sh 17 ./paralelo/fase2/bin card 2 >> resultadosFase2Jon.log
done
echo "> T 65 procesos grupos de 4" >> resultadosFase2Jon.log
for a in {1..3}; do
   exec.sh 65 ./paralelo/fase2/bin card 16 >> resultadosFase2Jon.log
done

echo "> T 65 procesos grupos de 8" >> resultadosFase2Jon.log
for a in {1..3}; do
   exec.sh 65 ./paralelo/fase2/bin card 8 >> resultadosFase2Jon.log
done

echo "> T 65 procesos grupos de 16" >> resultadosFase2Jon.log
for a in {1..3}; do
   exec.sh 65 ./paralelo/fase2/bin card 4 >> resultadosFase2Jon.log
done
echo "> T 112 procesos grupos de 16" >> resultadosFase2Jon.log
for a in {1..3}; do
   exec.sh 113 ./paralelo/fase2/bin card 7 >> resultadosFase2Jon.log
done
echo "> T 112 procesos grupos de 28" >> resultadosFase2Jon.log
for a in {1..3}; do
   exec.sh 113 ./paralelo/fase2/bin card 4 >> resultadosFase2Jon.log
done
echo "> T Resultados Fase2 Ivan" > resultadosFase2Ivan.log
echo "> T 17 procesos grupos de 4" >> resultadosFase2Ivan.log
for a in {1..3}; do
   exec.sh 17 ./paraleloIvan/bin card 4 >> resultadosFase2Ivan.log
done
echo "> T 17 procesos grupos de 8" >> resultadosFase2Ivan.log
for a in {1..3}; do
   exec.sh 17 ./paraleloIvan/bin card 8 >> resultadosFase2Ivan.log
done
echo "> T 65 procesos grupos de 4" >> resultadosFase2Ivan.log
for a in {1..3}; do
   exec.sh 65 ./paraleloIvan/bin card 4 >> resultadosFase2Ivan.log
done

echo "> T 65 procesos grupos de 8" >> resultadosFase2Ivan.log
for a in {1..3}; do
   exec.sh 65 ./paraleloIvan/bin card 8 >> resultadosFase2Ivan.log
done

echo "> T 65 procesos grupos de 16" >> resultadosFase2Ivan.log
for a in {1..3}; do
   exec.sh 65 ./paraleloIvan/bin card 16 >> resultadosFase2Ivan.log
done
echo "> T 112 procesos grupos de 16" >> resultadosFase2Ivan.log
for a in {1..3}; do
   exec.sh 113 ./paraleloIvan/bin card 16 >> resultadosFase2Ivan.log
done
echo "> T 112 procesos grupos de 28" >> resultadosFase2Ivan.log
for a in {1..3}; do
   exec.sh 113 ./paraleloIvan/bin card 28 >> resultadosFase2Ivan.log
done
