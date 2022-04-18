# Glosario de terminos


- Infiltración: Agua que infiltra desde las subcuencas de SWMM, hacia el aquífero (mundo subterraneo). Esto desaparece del dominio de SWMM para efectos de su simulación, para nosotros, esto se convierte en RECARGA (asociado al paquete `.rch`) en el modelo MODFLOW.
- Seepage: Término usado por SWMM para referirse a la infiltración que se da en las unidades de almacenamiento de SWMM (storage units)
- Escorrentía: Exceso de agua que existe en una subcuenca por lo que fluye (escurre/drena) a un nodo superficial (storage units, junctions, dividers, outfalls)
- Exfiltración: Agua que drena desde las aguas subterraneas hacia la superficie. Se calcula mediante el paquete `.drn` de MODFLOW y el volumen total drenado se incorpora en nodos superficiales como un flujo lateral (lateral inflow para SWMM).
- Lateral inflow: Parámetro de los nodos (storage units, junctions, dividers, outfalls) de SWMM. Este indica una tasa de flujo que ingresa al nodo. Nosotros lo usaremos para incorporar el agua exfiltrada, EN CADA PASO DE TIEMPO.