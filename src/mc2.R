#### Cleaning the environment ####
rm(list = ls())
gc()

#### Libraries ####
set.seed(123)
source('functions/loadlib.R')
libraries = c('tidyverse', 'lubridate', 'robust', 'xts', 'mgcv', 'nlme', 'RColorBrewer', 'plotly', "cluster")
for (i in libraries) loadlib(i)
rm('i', 'libraries')

#### Functions ####
procesar_clusters = function(dataset, clus_opt = "kmeans", rango_clus = 1:10, iter_max = 50) {
  if(clus_opt == "kmeans") {
    silh_list_ext = list()
    silh_avgwidth_ext = array()
    for (i in rango_clus) {
      print(i) #esto es para ir monitoreando el avance cuando corre
      modelo = kmeans(dataset, centers = i+1, algorithm = "MacQueen",
                          iter.max = iter_max) #calculo el kmeans
      # i          =  i + 10*(j-1) #esto es una correccion que le hago al i para que no se me sobreescriban en la lista. Esta funcion en j=2 vale 10, en j=3 vale 20 y as?
      distancias = dist(dataset, "euclidean")
      silh_list_ext[[i]] = silhouette(modelo$cluster, distancias) #guarda en una lista el silhouette
      silh_avgwidth_ext[i]   = summary(silh_list_ext[[i]])[["avg.width"]] #guarda el silhouette promedio
    }
  }
  else {
    stop("Por ahora solo kmeans es el unico algoritmo implementado")
    return(-1)
  }
  return(data_frame(silhouette_promedio = silh_avgwidth_ext,
                      cant_clusters = rango_clus+1))
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
  mutate(mes  = str_pad(month(`sample_date`), width = 2, pad = 0, side = "left"),
         year = year(`sample_date`)) %>%
  unite(col = "ym", year, mes, sep = "")

# Medidas
sort(unique(data.readings[["measure"]]))
sort(unique(data.readings[["measure"]])) %>% length(.) # 106 variables en total

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
rlinear.model = robust::lmRob(data = data.readings[data.readings[["measure"]] == "Methylosmoline", ], value ~ `sample_date`)

# TODO: tomar la media y la mediana por mes entre las variables en total y por estacion, para ver cambios en el parque.
data.month = data.readings %>%
    mutate(month_date = str_pad(string = month(sample_date), width = 2, pad = '0', side = 'left'),
           monthyear  = paste0(year(sample_date), month_date)) %>%
    group_by(monthyear, location, measure) %>%
    summarise(avg_val = mean(value, na.rm = TRUE),
              median_val = median(value, na.rm = TRUE),
              std_err = sd(value, na.rm = TRUE))

# TODO: lo que habría que hacer es clusterizar por variable y por estacion. Solamente asi se pueden detectar outliers entre las estaciones
# Por ejemplo, clusterizando para Methylosmoline
subset.methyl.kohsoom = data.month[data.month[["location"]] == "Kohsoom" & data.month[["measure"]] == "Methylosmoline", "median_val"]
kmeans1 = procesar_clusters(dataset = subset.methyl.kohsoom)

# Alrededor de 4 son los clusters que me da con silhouette alto
dist.methyl.kohsoom  = dist(subset.methyl.kohsoom, method = "euclidean")
model.methyl.kohsoom = kmeans(subset.methyl.kohsoom, centers = 4, iter.max = 50, algorithm = "MacQueen")

# Silhouette de este metodo
silhouette(x = model.methyl.kohsoom$cluster, dist.methyl.kohsoom)

data.year = data.readings %>%
  mutate(year_date = str_pad(string = year(sample_date), width = 2, pad = '0', side = 'left')) %>%
  group_by(year_date, location, measure) %>%
  summarise(avg_val = mean(value, na.rm = TRUE),
            median_val = median(value, na.rm = TRUE),
            std_err = sd(value, na.rm = TRUE))

# Clustering sobre los meses. Aca, asi como esta, puedo ver como se agrupan las distintas estaciones en cuanto a su mediana
# Pruebo hacer un cluster por meses por estacion, usando los quimicos como variable
data.month.wide = data.month %>%
    select(measure, location, monthyear, median_val) %>%
    spread(key = measure, value = median_val)

data.year.wide = data.year %>%
    select(measure, location, year_date, median_val) %>%
    spread(key = measure, value = median_val)

# K-means da error porque hay muchos huecos entre variables medidas en el tiempo, y huecos temporales al interior de cada variable también

# TODO: ¿Se puede detectar anomalias sin tener en cuenta la variable del tiempo?

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
