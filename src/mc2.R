#### Cleaning the environment ####
rm(list = ls())
gc()

#### Libraries ####
set.seed(123)
source('functions/loadlib.R')
libraries = c('tidyverse', 'lubridate', 'robust', 'xts', 'mgcv', 'nlme', 'RColorBrewer', 'plotly', 'cluster', 'factoextra')
for (i in libraries) loadlib(i)
rm('i', 'libraries')

#### Functions ####
procesar_clusters = function(dataset, clus_opt = "kmeans", rango_clus = 1:10, iter_max = 50) {
  if (clus_opt == "kmeans") {
    silh_list_ext = list()
    silh_avgwidth_ext = array()
    for (i in rango_clus) {
      print(i) #esto es para ir monitoreando el avance cuando corre
      modelo = kmeans(dataset, centers = i + 1, algorithm = "MacQueen",
                          iter.max = iter_max) #calculo el kmeans
      # i          =  i + 10*(j-1) #esto es una correccion que le hago al i para que no se me sobreescriban en la lista. Esta funcion en j=2 vale 10, en j=3 vale 20 y as?
      distancias = dist(dataset, "euclidean")
      silh_list_ext[[i]] = silhouette(modelo$cluster, distancias) #guarda en una lista el silhouette
      # fviz_silhouette(silh_list_ext[[i]], main = paste0("Silhouette de ", i, " clusters"))
      plot(silh_list_ext[[i]], color = "blue")
      silh_avgwidth_ext[i]   = summary(silh_list_ext[[i]])[["avg.width"]] #guarda el silhouette promedio
    }
  }
  else {
    stop("Por ahora solo kmeans es el unico algoritmo implementado")
    return(-1)
  }
  return(data_frame(silhouette_promedio = silh_avgwidth_ext,
                      cant_clusters = rango_clus + 1))
}

#### Variables ####
readings = "data/Boonsong Lekagul waterways readings.csv"
units    = "data/chemical units of measure.csv"

#### Datasets ####
data.readings = read_csv(readings, locale = locale(encoding = "latin1"))
data.units    = read_csv(units, locale = locale(encoding = "latin1"))

data.readings = rename(data.readings, "sample_date" = "sample date")

#### Program ####
# uniendo con las unidades de medidas
data.readings = data.readings %>%
  left_join(data.units, by = "measure")

data.readings[["sample_date"]] = dmy(data.readings[["sample_date"]])

# Creando una variable de periodo: monthyear
data.readings = data.readings %>%
  mutate(mes  = str_pad(month(sample_date), width = 2, pad = 0, side = "left"),
         year = year(sample_date)) %>%
  unite(col = "ym", year, mes, sep = "")

# Medidas
sort(unique(data.readings[["measure"]]))
sort(unique(data.readings[["measure"]])) %>% length() # 106 variables en total

# Cantidad de muestras por mes de methylosmoline por estación
# El methylosmoline es, supuestamente, el químico que mayormente afecta al pipit.
data.readings %>%
  filter(measure == "Methylosmoline") %>%
  group_by(location, ym) %>%
  summarise(value = sum(value),
            cuenta = n()) %>%
  ggplot(data = ., aes(x = ym, y = location, fill = cuenta)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cantidad de muestras de methylosmoline tomadas por mes", x = "AñoMes", y = "Estación", fill = "Cantidad") +
  geom_tile()

# Boxplot of measures
boxplot(data.readings$value[data.readings$measure == "Methylosmoline"] ~ data.readings$location[data.readings$measure == "Methylosmoline"])

# Methylosmoline
data.readings %>%
  filter(measure == "Methylosmoline") %>%
  group_by(`sample_date`, location, measure) %>%
  transmute(value = mean(value)) %>%
  ggplot(data = ., aes(x = `sample_date`, y = value, color = location)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 150, by = 15)) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  theme_bw() +
  labs(title = "Methylosmoline", x = "Fecha", y = "Promedio de valores de concentración (µg/l)")

linear.model  = lm(data = data.readings[data.readings[["measure"]] == "Methylosmoline", ], value ~ `sample_date`)
rlinear.model = robust::lmRob(data = data.readings[data.readings[["measure"]] == "Methylosmoline", ], value ~ sample_date)

data.month = data.readings %>%
    mutate(month_date = str_pad(string = month(sample_date), width = 2, pad = '0', side = 'left'),
           monthyear  = paste0(year(sample_date), month_date)) %>%
    group_by(monthyear, location, measure) %>%
    summarise(avg_val = mean(value, na.rm = TRUE),
              median_val = median(value, na.rm = TRUE),
              std_err = sd(value, na.rm = TRUE))

# Tomo las variables de un ejemplo: una de las estaciones
data.month.kohsoom = data.month %>%
  filter(location == "Kohsoom") %>%
  select(measure, median_val) %>%
  spread(key = measure, value = median_val)

data.month.kohsoom

# TODO: lo que habría que hacer es clusterizar por variable y por estacion. Solamente asi se pueden detectar outliers entre las estaciones
# Por ejemplo, clusterizando para Methylosmoline
subset.methyl.kohsoom = data.month[data.month[["location"]] == "Kohsoom" & data.month[["measure"]] == "Methylosmoline", "median_val"]
kmeans1 = procesar_clusters(dataset = subset.methyl.kohsoom)

# Alrededor de 4 son los clusters que me da con silhouette alto
dist.methyl.kohsoom  = dist(subset.methyl.kohsoom, method = "euclidean")
model.methyl.kohsoom = kmeans(subset.methyl.kohsoom, centers = 4, iter.max = 50, algorithm = "MacQueen")

# Silhouette de este metodo
silhouette(x = model.methyl.kohsoom$cluster, dist.methyl.kohsoom)

# Prueba el cluster que supongo que hizo Florencia
kmeans2 = procesar_clusters(dataset = data.month[data.month[['measure']] == 'Total dissolved salts', 'median_val'], rango_clus = 1:14)

data.year = data.readings %>%
  mutate(year_date = str_pad(string = year(sample_date), width = 2, pad = '0', side = 'left')) %>%
  group_by(year_date, location, measure) %>%
  summarise(avg_val = mean(value, na.rm = TRUE),
            median_val = median(value, na.rm = TRUE),
            std_err = sd(value, na.rm = TRUE))

# Prueba el cluster que supongo que hizo Florencia
kmeans3 = procesar_clusters(dataset = data.year[, 'median_val'], rango_clus = 1:14)

# Parece que no es exactamente lo que estoy probando, habria que modificarlo.
# TODO: hay que armar una prediccion de valores en base a la historia. Tiene que evitar problemas como el tema de que se agreguen ceros
# TODO: tambien pueden armarse clusters entre variables en una locacion, o entre variablse de distintas locaciones de, por ejemplo, una misma cuenca.
# TODO: Como herramienta, estaria bueno que puedan manipularse las estaciones por cuenca.
# TODO: Lo que voy a intentar hacer es armar una prediccion de outliers basado en el rango de las variables
# TODO: tendria que aceptar tambien el input del investigador, por medicion, por lugar y por momento en el tiempo, ya que los rangos aceptados deberian poder cambiar

# Clustering sobre los meses. Aca, asi como esta, puedo ver como se agrupan las distintas estaciones en cuanto a su mediana
# Pruebo hacer un cluster por meses por estacion, usando los quimicos como variable
data.month.wide = data.month %>%
    select(measure, location, monthyear, median_val) %>%
    spread(key = measure, value = median_val)

data.year.wide = data.year %>%
    select(measure, location, year_date, median_val) %>%
    spread(key = measure, value = median_val)

# K-means da error porque hay muchos huecos entre variables medidas en el tiempo, y huecos temporales al interior de cada variable también

# Clusterizacion a medida del usuario
year_inicio   = '2014'
mes_inicio    = '01'
loc           = 'Kohsoom'

data.kmeans = data.month.wide[data.month.wide[["location"]] == loc,]

# Cantidad de nulos por variable
filtro_na = sapply(data.kmeans, function(x) sum(is.na(x))) %>% sort
na.permitidos = 10

# Variables elegidas segun la cantidad de nulos
variables.elegidas = names(filtro_na)[filtro_na <= na.permitidos]
variables.elegidas = variables.elegidas[!(variables.elegidas %in% c("location", "monthyear"))]

# TODO: en este momento es que el usuario deberia poder ver las variables que quedaron tras el filtro de nulos. Aca debería poder elegir a mano con las que quisiera quedarse o sacar.

# Filtrando el dataset segun las variables elegidas
data.kmeans = data.kmeans[, c("location", "monthyear", variables.elegidas)]

# Por default, la herramienta hace un escalamiento común.
robust.scale = FALSE

# Escalando las variables segun su rango. Puede ser un escalamiento robusto (restando mediana y dividiendo por MAD) o escalamiento estandar (restando media y dividiendo por desvio)
if (robust.scale) {
  loadlib("quantable")
  data.kmeans[, variables.elegidas] = robustscale(data.kmeans[, variables.elegidas])$data
} else {
  for (i in variables.elegidas) data.kmeans[[i]] = scale(data.kmeans[[i]])
}

# TODO: dejarle al usuario elegir la cantidad de clusters a probar. Se recomienda por lo menos de 2 a 10 clusters como default.

# Clusterización: el usuario puede hacerla a su medida.
# Por default, se prueba de 2 a 10 clusters
cant.clusters.prueba = 1:9

# TODO: de los clusters que elija, mostrar gráficos de las metricas de silhouette promedio, las sumas de cuadrados entre clusters y dentro de cada cluster.
# TODO: elegir el cluster segun la mejor metrica de silhouette promedio. No obstante, dejar al usuario cambiar la cantidad de clusters que quiere calcular, mostrando un gráfico de silhouettes por cluster para ayudarlo a decidir: lo mejor sería que pueda comparar para ver si elige alguna agrupacion que le permita un cluster bien discriminado con alto silhouette a pesar de que el silhouette promedio no sea el mejor de todos.
# TODO: regresiones: implementar las regresiones normales, polinómicas, smoothing y regresiones robustas (ver si hay robustas más que lineales) para detección de outliers. Implementar detección de outliers de acuerdo a lo que vimos con Soria.

# Resolución del problema del VAST challenge.
# TODO: ¿Qué más se puede aportar respecto a lo que presentamos en el challenge, con relación al Methylosmoline y a cómo afecta al Pipit? Rankear las estaciones por cuenca y según la distancia a la supuesta zona de tirada de desechos.
# TODO: Desafíos a futuro para la herramienta: lo ideal sería que se pueda seguir usando la herramienta incorporando técnicas de series de tiempo y correlación espacio-temporal entre las estaciones, por compartir la misma cuenca, el mismo parque y por la distancia a la zona donde se arrojan los desechos.
# TODO: ¿hay algún cambio distribucional en el tiempo de los desechos? ¿Y con respecto al espacio? ¿hay alguna expansión de la marea de químicos (posiblemente el methylosmoline o altamente correlacionados con este químico) en el espacio lejos de la zona de desechos, en el tiempo?


# Poniendo las variables como columnas, en vez de key-value
spread.data = spread(data = data.readings, key = c("sample_date", "measure"), value = "value")

# Poniendo las variables como columnas, en vez de key-value
spread.data = spread(data = data.readings, key = c("sample_date", "measure"), value = "value")
# Frecuencia de muestreo de las variables
names(data.readings)[4] = "date"

spread.vars = data.readings %>%
  group_by(date, measure, location) %>%
  summarise(cuenta = n(),
            media = mean(value),
            mediana = median(value),
            maximo = max(value),
            minimo = min(value))

spread.vars = data.readings %>%
  spread(measure, value)

table(spread.vars$cuenta)

# Tratando de normalizar todas las variables entre 0 y 1.
split.data.readings = data.readings %>%
  group_by(measure) %>%
  mutate(normalizedValue = (value - min(value)) / max(value))

split.data.readings[["normalizedValue"]][is.na(split.data.readings[["normalizedValue"]])] = 0 # El berilio da siempre cero. Lo imputo a cero, por ende

# Habiendo normalizado las variables, quiero fijarme si hay un aumento generalizado del promedio de mediciones en algún momento del tiempo por estación.

split.data.readings %>%
  group_by(location, date) %>%
  summarise(valores = mean(normalizedValue)) %>%
  ggplot(aes(x = date, y = valores, colour = location)) +
  geom_line() +
  labs(title = "Mediciones en el tiempo por estación", x = "Fecha", y = "Valores")

# Separando las estaciones para mayor claridad
split.data.readings %>%
  group_by(location, date) %>%
  summarise(valores = mean(normalizedValue)) %>%
  ggplot(aes(x = date, y = valores)) +
  geom_line() +
  theme_bw() +
  labs(title = "Mediciones en el tiempo por estación", x = "Fecha", y = "Valores") +
  facet_wrap(~location)

# Lo primero que se observa en este gráfico es que no todas las estaciones tuvieron observaciones en el mismo periodo de tiempo, por lo que no son tan comparables las series de tiempo, sobre todo si queremos obtener correlaciones espacio-temporales entre estaciones.

# Medidas por año
graph = data.readings %>%
  filter(measure == "Water temperature") %>%
  ggplot(aes(x = as.character(year(date)), y = value)) +
  theme_bw() +
  geom_boxplot() +
  labs(title = "Temperatura del agua", x = "Año", y = "Valores") +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)
graph = ggplotly(graph)
graph

# "Iron" parece tener un problema en el año 2003.
graph = data.readings %>%
  filter(measure == "Iron") %>%
  ggplot(aes(x = as.character(year(date)), y = value)) +
  theme_bw() +
  geom_boxplot() +
  labs(title = "Hierro", x = "Año", y = "Valores") +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)
graph = ggplotly(graph)
graph

# Si saco ese año del ploteo
graph = data.readings %>%
  filter(measure == "Iron" & year(date) != 2003) %>%
  ggplot(aes(x = as.character(year(date)), y = value)) +
  theme_bw() +
  geom_boxplot() +
  labs(title = "Hierro", x = "Año", y = "Valores") +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)
graph = ggplotly(graph)
graph

# Sacando los outliers más pronunciados, se detectan otros outliers distintos en el rango de las variables.
# ¿Puede ser que Iron tenga, en el 2003, un problema sencillo de puntos y comas (separadores decimales y de miles)?

graph = data.readings %>%
  filter(measure == "Methylosmoline" & year(date) != 2003) %>%
  ggplot(aes(x = as.character(year(date)), y = value)) +
  theme_bw() +
  geom_boxplot() +
  labs(title = "Methylosmoline", x = "Año", y = "Valores") +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)
graph = ggplotly(graph)
graph

rm(spread.vars2)
gc()

#### Clusterizacion ####
# Lo que voy a hacer es armar las variables agregadas por mes con las mediciones como variables.
# La agregacion que voy a hacer por mes es la mediana de las mediciones.
# La idea es ver si varias de las medidas se mueven juntas.
