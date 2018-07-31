#### Cleaning the environment ####
rm(list = ls())
gc()

#### Libraries ####
source('functions/loadlib.R')
libraries = c('tidyverse', 'lubridate', 'robust', 'xts', 'mgcv', 'nlme', 'RColorBrewer', 'plotly')
for (i in libraries) loadlib(i)
rm('i', 'libraries')

#### Variables ####
readings = "data/Boonsong Lekagul waterways readings.csv"
units    = "data/chemical units of measure.csv"

#### Datasets ####
data.readings = read_csv(readings, locale = locale(encoding = "latin1"))
data.units    = read_csv(units, locale = locale(encoding = "latin1"))

#### Program ####
# Uniendo con las unidades de medidas
data.readings = data.readings %>% 
  left_join(data.units, by = "measure")

data.readings[["sample date"]] = dmy(data.readings[["sample date"]])

# Creando una variable de periodo: monthyear
data.readings = data.readings %>% 
  mutate(mes  = str_pad(month(`sample date`), width = 2, pad = 0, side = "left"),
         year = year(`sample date`)) %>% 
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
  group_by(`sample date`, location, measure) %>% 
  transmute(value = mean(value)) %>% 
  ggplot(data = ., aes(x = `sample date`, y = value, color = location)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 150, by = 15)) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  theme_bw() +
  labs(title = "Methylosmoline", x = "Fecha", y = "Promedio de valores de concentración (µg/l)")

data.readings %>% 
  filter(measure == "p,p-DDT") %>% 
  group_by(`sample date`, location, measure) %>% 
  transmute(value = mean(value)) %>% 
  ggplot(data = ., aes(x = `sample date`, y = value, color = location)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 2, by = 0.2)) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  labs(title = "p,p-DDT", x = "Fecha", y = "Promedio de valores de concentración (µg/l)")

data.readings %>% 
  filter(measure == "p,p-DDE") %>% 
  ggplot(data = ., aes(x = `sample date`, y = value, color = location)) +
  geom_point() +
  geom_line()

linear.model  = lm(data = data.readings[data.readings[["measure"]] == "p,p-DDT", ], value ~ `sample date`)
rlinear.model = robust::lmRob(data = data.readings[data.readings[["measure"]] == "p,p-DDT", ], value ~ `sample date`)

plot(x = data.readings[['']][data.readings[['measure']] == 'p,p-DDT'],
     y = data.readings[['value']][data.readings[['measure']] == 'p,p-DDT'])

spread.data = spread(data = data.readings, key = c("sample date", "measure"), value = "value") 

# Frecuencia de muestreo de las variables
names(data.readings)[4] = "date"

spread.vars = data.readings %>% 
  group_by(date, measure, location) %>%
  summarise(cuenta = n(),
            media = mean(value),
            mediana = median(value),
            maximo = max(value),
            minimo = min(value))

spread.vars2 = data.readings %>%
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
