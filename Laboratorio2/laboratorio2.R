library(tidyverse)
data(starwars)

#select →Mantener o eliminar columnas utilizando sus nombres y tipos
#filter → Mantener filas que coincidan con una condición
#group_by → Agrupar por una o más variables
#tally →Cuenta las observaciones en cada grupo
#summarise → Resume cada grupo en una fila y puede usarse para calcular parámetros
#estadísticos con funciones básicas [mean() ; median (); sd() ; min() ; max() ; IQR()]
#read_csv → Leer un archivo delimitado (incluidos CSV y TSV) en un data frame



#Ejemplos usando select
  # Seleccionar todas las columnas menos el nombre
  starwars %>% select(-name)
  #Seleccionar sólo las columnas que tienen subraya (_)
  starwars %>% select(contains("_"))
  #Seleccionar sólo las columnas que empiezan con "s"
  starwars %>% select(starts_with("s"))
  #Crear un data frame con los nombres y planeta de origen (homeworld)
  homeworld <- starwars %>% select(name, homeworld)
  #Filtrar datos 
    #Filtrar por especies: sólo humanos
    human <- starwars %>% filter(species == "Human")
    #Filtrar por especies: sólo humanos del planeta Tatooine
    starwars %>% filter(species == "Human", homeworld == "Tatooine")
    #Crear un nuevo datframe con todas las especies menos los Droides
    starwars_nodroids <- starwars %>% filter(species != "Droid")
      #¿Cuántos registros cumplen las condiciones finales?
        #En homeworld 87
        #En human 35
        #En starwars 87 
        #En starwars_nodroids 77

#Seleccionar y agrupar datos
  #Usamos group_by y tally
  starwars %>% group_by(species) %>% tally()
  #Añadiendo otra variable
  starwars %>% group_by(species, gender) %>% tally()
  #Si lo quieres guardar en el environment recuerda asignarle un nombre
  table_gender <- starwars %>% group_by(species, gender) %>% tally()

#Calcular algunos estadísticos MEDIA Y DESVIACION ESTANDAR
  #a.rm=T quiere decir que elima los NA (valores No Asignados o sin datos) Asi calculamos la media
  starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm = T),mean_mass = mean(mass,na.rm = T))
  #¿Cómo calcularías la desviación estándar (sd) de esos parámetros? 
    #group_by(species) nos lo separa por especies
    starwars %>% group_by(species) %>% summarise(sd_height = sd(height, na.rm = T),sd_mass = sd(mass,na.rm = T))
    
#Crear gráficos y modificar algunos elementos.
  #Hacer un gráfico de la altura vs. la masa de los personajes
  ggplot(starwars, aes(height, mass)) + geom_point()
  #Puedes modificar el color 
  ggplot(starwars, aes(height, mass)) + geom_point(colour = "red")
  #Modificando el color y el punto
  ggplot(starwars, aes(height, mass)) + geom_point(colour = "purple", pch = 3)
  #Modificando el color y el fondo 
  ggplot(starwars, aes(height, mass)) + geom_point(colour = "red") + theme_light()
  #Al crear los gráficos puedes observar que hay un punto que corresponde a un personaje con una masa muy grande. Inspecciona el datset, filtra usando las funciones de tidyverse, crea un nuevo dataframe sin ese personaje y crea de nuevo el gráfico final. (Exporta el gráfico con la opción exportar en el panel derecho y adjúntalo en el pdf)
    #Podemos filtrar los valores de mass que sean anormalmente altos para encontrar al personaje problemático. 
    # Inspeccionar los valores altos de masa
    starwars %>% filter(mass > 500)  # Filtra personajes con masa mayor a 500. En este caso es Jabba
    #Lo eliminamos del dataset usando filter(), y creamos un nuevo dataframe sin ese personaje. Filtrar el personaje con masa anómala (mayor de 500)
    starwars_clean <- starwars %>% filter(mass < 500)
    # Crear el gráfico sin el personaje con masa anómala
    ggplot(starwars_clean, aes(height, mass)) + 
      geom_point(colour = "red") + 
      theme_light()
    
#Descarga el dataset toy.csv cargalo en R studio usando la función read_csv de la libreria tidyverse. Tienes que poner la dirección donde has guardado el archivo descargado. En el ejemplo, el archivo está en la carpeta "Descargas"
  #ir a descargas y copiar ruta de acceso, si no va cambiar barras
  #Si no va dar a import data base
  toy <- read_csv("C:/Users/ruizl/Downloads/toy.csv")
  #IMC es Indice de Masa Corporal, IMC_clas: clasificación de acuerdo al IMC.
  #IAS significa Indice de Alimentación Saludable; IAS_clas: clasificación de acuerdo al IAS.
  #CCintura significa Circunferencia de cintura.
      #Inspecciona el dataset, haz un resumen de la media (mean) de las variables (Peso, Altura,IMC, IAS, CCintura). Agrupando por sexo.
          #1. INSPECCION Ver las primeras filas del dataset
          head(toy)
            # Ver la estructura del dataset
            str(toy)
            # Resumen estadístico de las variables
            summary(toy)
          #2. Resumen de la media variables
            # Calcular la media de las variables agrupadas por sexo
            toy_summary <- toy %>% group_by(Sex) %>% summarise(mean_peso = mean(Weight_Kg, na.rm = TRUE),mean_altura = mean(Height_cm, na.rm = TRUE),mean(IMC, na.rm = TRUE),mean_IAS = mean(IAS, na.rm = TRUE), mean_CCintura = mean(Ccintura, na.rm = TRUE)) #Poner exactamente como lo pone en la tabla
            # Mostrar el resumen
            toy_summary
         #3. Haz una tabla sólo con los pacientes femeninos ¿Cuántos registros cumplen las condiciones? ¿De estos cuantos tienen Sobrepeso (Overweight)?  Usa select y filter.
            # Filtrar pacientes femeninos
            pacientes_femeninos <- toy %>% filter(Sex == "Women") %>% select(Sex, Age, Weight_Kg, Height_cm, IMC, IMC_clas, IAS, Ccintura, IAS_Clas)
            pacientes_femeninos
            # Contar cuántos registros son de pacientes femeninos
            num_femeninos <- nrow(pacientes_femeninos) #58
            # Filtrar pacientes femeninos con sobrepeso (Overweight)
            pacientes_sobrepeso <- pacientes_femeninos %>% filter(IMC_clas == "Overweight")
            # Contar cuántos registros tienen sobrepeso
            num_sobrepeso <- nrow(pacientes_sobrepeso) #9
            # Mostrar resultados
            num_femeninos   # Número total de pacientes femeninos
            num_sobrepeso   # Número de pacientes femeninos con sobrepeso
        #4.  Haz un gráfico usando ggplot relacionando el IMC (Indice de masa corporal) con el peso (Weight_Kg) de todos los pacientes.
             ggplot(toy, aes(x = Weight_Kg, y = IMC)) + geom_point(color = "blue") + labs(title = "Relación entre IMC y Peso", x = "Peso (Kg)", y = "IMC (Índice de Masa Corporal)")
              #ggplot(toy, aes(x = Weight_Kg, y = IMC)): Esto inicializa el gráfico usando el dataset toy, definiendo Weight_Kg como el eje x y IMC como el eje y.
              #geom_point(color = "blue"): Añade un gráfico de dispersión (scatter plot) con los puntos en color azul.
              #labs(): Añade etiquetas personalizadas al gráfico, como el título y los nombres de los ejes.
        #5.Repítelo filtrando sólo los pacientes categorizados como "Overweight" y "Obesity".
            #filter(IMC_clas %in% c("Overweight", "Obesity")): Filtra los datos para incluir solo a los pacientes que están clasificados como Overweight o Obesity en la columna IMC_clas.
            ggplot(toy %>% filter(IMC_clas %in% c("Overweight", "Obesity")), aes(x = Weight_Kg, y = IMC)) + geom_point(color = "blue") + labs(title = "Relación entre IMC y Peso (Overweight y Obesity)", x = "Peso (Kg)", y = "IMC (Índice de Masa Corporal)")

#Ejercicio:Utiliza los comandos adecuados para instalar los paquetes de R ape phangorn y phytools que utilizaremos en el laboratorio de la siguiente semana. Carga las librerías y envia un print de pantalla con el output, demostrando que la instalación fue exitosa o si hubo algún problema.
            install.packages("ape")
            install.packages("phangorn")
            install.packages("phytools")
            