# Explorando la independencia espacial en datos de cámaras trampa
# Gabriel Andrade Ponce


# Librerias ---------------------------------------------------------------

library(leaflet)# Create Interactive Web Maps with the JavaScript 'Leaflet'
library(sp)# Classes and Methods for Spatial Data
library(tidyverse)# Easily Install and Load the 'Tidyverse'
library(AICcmodavg) # Model Selection and Multimodel Inference Based on (Q)AIC(c)
library(DHARMa) # Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models
library(ncf) # Spatial Covariance Functions
library(MASS) # Support Functions and Datasets for Venables and Ripley's MASS
library(showtext) # Using Fonts More Easily in R Graphs
font_add_google("Special Elite", family = "special") # Tipo de letra
showtext_auto() #llamar la letra


# Mapa --------------------------------------------------------------------

CToperation <- read.csv("Data/CTtable.csv") 

# Generar objeto espacial para el mapa
CT_points <- SpatialPoints(cbind(CToperation$utm_x, 
                                 CToperation$utm_y),
                           proj4string = CRS('+proj=utm +datum=WGS84 +zone=14 +towgs84=0,0,0')) 

# Proyectar a WGS solo para este paso
CT_points <- spTransform(CT_points, "+proj=longlat +datum=WGS84")

# Generar el mapa
m <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group="Satellite") %>%  # Add satellite data
  addProviderTiles(providers$Esri.WorldTopoMap, group="Base") %>% 
  addCircleMarkers(lng=coordinates(CT_points)[,1], lat=coordinates(CT_points)[,2], 
                   popup= paste(CToperation$Station)) %>% 
  # Layers control
  addLayersControl(
    baseGroups = c("Satellite", "Base"),
    options = layersControlOptions(collapsed = FALSE)
  )
m


# Cargar datos ------------------------------------------------------------

# Tabla de funcionamiento de cámaras trampa con coordenadas
CToperation <- read.csv("Data/CTtable.csv") 

# Usaremos el archivo que creamos en survey report, se llama events_by_station2.csv

freq_reg <- read.csv("Data/surveyReport/events_by_station2.csv")%>% # Llamamos el csv
  filter(Species== "Odocoileus virginianus") %>%  # Usaremos solo los datos de venados
  left_join(CToperation, by= "Station") # Unimos para agregar las coordenadas


# Mapa de registros -------------------------------------------------------

(mapa <- freq_reg %>% # Llamamos los datos
    ggplot(aes(x=utm_x, y=utm_y))+ # Establecemos x y y
    geom_point(aes(size= n_events), alpha=0.9, colour= "steelblue")+ # Decimos que queremos geom de puntos
    scale_size(range = c(1,15))+ # La escala de los puntos
    labs(title= "Mapa de número detecciones de venado cola blanca", 
         size= "Número de \ndetecciones")+
    theme_bw()+ # Tema
    theme(text=element_text(size=15, family = "special"))) #tipo de letra


# Modelos lineales generalizados para la frecuencia de registro -----------
# Generamos el data frame
sp_glmdata <- read.csv("Data/covs.csv") %>% # Llamamos las covariables que vamos a usar
  right_join(freq_reg, by= "Station")# Las unimos con nuestra tabla de número de eventos

# Modelos lineales generalizados 

# sin variables
m0 <- glm(n_events~ 1, # Formula
          data = sp_glmdata, # datos
          family = "poisson") # tipo de distribución

# la frecuencia de registro afectada por la distancia a cultivo
m1 <- glm(n_events~ Dcrops, 
          data = sp_glmdata, 
          family = "poisson")

# la frecuencia de registro afectada por el verdor de la vegetación
m2 <- glm(n_events~ MSAVI, 
          data = sp_glmdata, 
          family = "poisson")

# la frecuencia de registro afectada por la pendiente
m3 <- glm(n_events~ Slope, 
          data = sp_glmdata, 
          family = "poisson") 

# la frecuencia de registro afectada por la distancia a poblados
m4 <- glm(n_events~ Dpop_G, 
          data = sp_glmdata, 
          family = "poisson")

# la frecuencia de registro afectada por el tipo de habitat
m5 <- glm(n_events~ Habitat, 
          data = sp_glmdata, 
          family = "poisson")

# para inspeccionar el modelo usar la función summary

summary(m5)

# Selección de modelos ----------------------------------------------------

lista_mods <- list(m0, m1, m2, m3, m4, m5) # crear lista de modelos
# Generar objeto con nombres de cada modelo
mod_names <- c("freq~ 1",
               "freq~ D_cultivos",
               "freq~ MSAVI",
               "freq~ Slope",
               "freq~ D_poblado",
               "freq~ Habitat"
)

AIC <- aictab(lista_mods, # Lista de modelos
              modnames = mod_names,  # nombres
              second.ord = F, # Seleccionar AIC
              sort = T) # Ordenar por menor valor

View(AIC)


# Verificando mejor modelo ------------------------------------------------

# Debido a que para glm poisson los residuales no se definen directamente, usamos simulateResiduals del paquete DHARMa
residuales <- simulateResiduals(fittedModel = m5, plot =F)

# Verificamos visualmente que el modelo cumpla los requisitos de la distribución
plotQQunif(residuales)


# Verificando autocorrelación espacial ---------------------------------------------

# Creamos el data.frame para el correlograma
data_resm5 <- data.frame(res=residuals(residuales), # Residuales que creamos
                         x= sp_glmdata$utm_x,# coordenadas en x
                         y= sp_glmdata$utm_y) # coordenadas en y


# Usamos la función correlog() del paquete nfc

m5_cor <- correlog(x=data_resm5$x, # coordenadas en x
                   y=data_resm5$y, # coordenadas en y
                   z=data_resm5$res, # variable de interés
                   na.rm=T, # en caso de NAs
                   increment = 400, # Distancia mínima de unidades
                   resamp=500) # número de iteraciones (se recomiendan mil)

# plot(m5_cor)
# puedes obtener el correlograma con plot(m5_cor), o hacer la gráfica en ggplot siguiendo el siguiente código

# Modificamos data frame para la gráfica
cor_m5 <- data.frame(correlation=m5_cor$correlation, 
                     distance= m5_cor$mean.of.class, 
                     p= m5_cor$p)%>% 
  mutate(p_valor= if_else(p<0.025, "significativo", "no-significativo"))

# Graficamos con ggplot

ggplot(cor_m5,aes(x=distance, y=correlation))+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_line( size=0.9)+
  geom_point(aes(colour= p_valor),size=3)+
  scale_x_continuous(breaks = seq(0,8000, by=500))+
  ylim(-1,1)+
  labs(x= "Unidades de distancia (m)", y= " Moran I", 
       title = " Correlograma residuales ")+
  theme_classic()+
  theme(text=element_text(size=15, family = "special" ),
        axis.text.x = element_text(angle = 45, hjust=1))



# Modelo con distribución binomial negativa -------------------------------
# Usando la función glm.nb de la paquetería MASS generamos un modelo con distribución binomial negativa

m5bn <- glm.nb(n_events~ Habitat, # Formula
               data = sp_glmdata) # datos

# Calculamos los residuales del nuevo modelo
residuales.bn <- simulateResiduals(fittedModel = m5bn, plot =F)

# Plot de los residuales
plotQQunif(residuales.bn)

# Fin --------------------------------------------------------