## Glossary (in Spanglish)

- **Infiltración:**
Dícese del agua que infiltra desde las subcuencas de SWMM hacia el [acuífero].
Para efectos de su simulación, esta desaparece del dominio de [SWMM] puesto que,
en el contexto de esta librería,
esta se convertirá en recarga asociada al paquete **[rch]** —de «recharge»— de [MODFLOW].

- **Exfiltración:**
Dícese del agua que drena desde las aguas subterráneas hacia la superficie.
Esta se calcula mediante el paquete **[drn]** —de «drain»— de MODFLOW.
El volumen total drenado se incorpora en nodos superficiales como flujo lateral (cf. *lateral inflow*) para SWMM.

- **Lateral inflow:**
Dícese de un parámetro de los nodos (e.g. *storage units*, *junctions*, *dividers*, *outfalls*)
de SWMM que indica una tasa de flujo que ingresa al nodo.
En esta librería se utiliza para incorporar el agua exfiltrada **en cada paso de tiempo**.

- **Escorrentía:**
Dícese del exceso de agua que existe en una subcuenca que termina drenando (o escurriendo) a un nodo superficial:
e.g. *storage units*, *junctions*, *dividers*, *outfalls*.

- **Seepage:**
Dícese de la infiltración que, en el contexto de SWMM, ocurre en las unidades de almacenamiento.

[modflow]:https://en.wikipedia.org/wiki/MODFLOW
[swmm]:https://en.wikipedia.org/wiki/Storm_Water_Management_Model

[drn]:https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/drn.html
[rch]:https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/rch.html

[acuífero]:https://es.wikipedia.org/wiki/Agua_subterránea
[escorrentía]:https://es.wikipedia.org/wiki/Escorrentía
