#LABORATORIO 3 - ARBOLES
  #1. Selección de los datos. (Tipo de secuencias: genes, proteínas…)
  #2. Alineamiento múltiple 
    # --> en vez de cuantificar el grado de similitud, se cuantifican las diferencias.
  #3. Definir modelo evolutivo que explique las diferencias observadas.
    #--> dan diferente valoracion a la recosntruccion para poder reconstruir mejor el arboles
  #4. Elegir método de construcción de árboles. Generar posibles árboles.ENRAIZAR
    #    --> Visualización y evaluación de los árboles.
  #5. Encontrar el mejor árbol y evaluar su fiabilidad. 

library(ape)
library(phangorn)
library(phytools)
#Primero nuevo proyecto con laborio tres y guardar en la carpeta donde este ese proyecto nuestas secuencias

#1. Alineamiento de la secuencia 
    #Descargar la secuencia 
    #Ve a ClustalX y realiza el alineamiento. Pon formato FASTA/PERSON.Darle a descargar.
    #Guardarlo en BLOC NOTAS llamado “fraxatin_aligned.fasta” (Recuerda guardalo en la carpeta donde está el proyecto de R)

#2.MATRIZ DISTANCIA - MODELO EVOLUTIVO
   #Para empezar convertimos el archivo en un formato compatible con la librería phangorn
   #Indicando el formato: “FASTA” y el tipo de secuencias: “AA” (aminoácidos)
   fraxatin <- read.phyDat(file = "fraxatin_aligned.fasta.txt", 
                         format = "FASTA", type = "AA")
  #Añadir el punto txt y si sale
  #Hay 11 secuencias con 310 caracteres y 180 columnas en las que hay diferentes posiciones. Los estados son los aminoácidos A (alanina), R (arginina), etc. 
  #Ahora lo que se hace es crear una matriz de distancia para poder crear árboles de distancia o parsimonia a través de ella. 
  #En este caso la clase debe ser AAbin (Amino Acid Sequences), por eso transformamos el objeto fraxatin en este tipo de clase.
      #La función dist.aa (método o función de objeto o clase) calcula una matriz de distancias por pares de las secuencias de aminoácidos partir de una objecto de clase AAbin utilizando un 
      #Utilizando modelo de evolución de aminoácidos (por ejemplo Dayhoff). Describe las probabilidades de que un aa en una secuencia cambia a otro aa durante el proceso de evolucion
      matrizdist <- as.AAbin(fraxatin)
      matrizdist <- dist.aa(matrizdist)
      matrizdist
      
#3.ARBOLES DE DISTANCIA - METODO CONSTRUCCION 
  #Árboles de similaridad, no representan inferencia evolutiva
    #METODO UPGMA 
        #Con la matriz de distancia perdemos los caracteres a favor de las diferencias en caracteres entre especies. 
        #Cuando aparece un valor de 0, significa que no hay diferencia al nivel de los caracteres (aminoácidos) entre las sequencias de dos especies. 
        #Ahora creamos un árbol con el método de grupo de pares no ponderados con media aritmética (UPGMA) usando la matriz de distancia que acabamos de calcular.
          arbolUPGMA <- upgma(matrizdist)
          plot(arbolUPGMA)
        #Si la longitud de dos ramas es indéntica, significa que las secuancias también son indénticas y en la matriz de distancia la diferencia es de 0.

    #METODO NJ -método de unión de vecinos (NJ) usando la misma matriz de distancias
          arbolNJ <- nj(matrizdist)
          plot(arbolNJ)
          #Este último árbol puede ser distinto del árbol creado con el método UPGMA. 
          #Para personalizar los árboles podemos agregar argumentos a parámetros como cex, para el tamaño de la letra, edge.color, para el grosos de las ramas.
          #También se puede escoger entre diferentes visualizaciones de árbol como filograma, cladograma, radial y demás.
          plot(arbolUPGMA, type= "p", cex=0.8, edge.width=2, edge.color="red", font=3)
          plot(arbolUPGMA, type= "c", cex=0.8, edge.width=2, edge.color="blue", font=3)
          plot(arbolUPGMA, type= "p", label.offset=0.0005, edge.lty=1, node.pos=2, cex=0.8, edge.width=2, edge.color="black", font=3)
              
          #PERSONALIZAR ARBOLES CON PHYTOOLS - GRAFICAR ARBOL
                #método plotTree del paquete phytools, el cual es compatible con ape y con phangorn.
                plotTree(arbolNJ)
                #A este también le podemos hacer modificaciones.
                plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="red", lwd=2)
                # En los árboles, sin cambiar la topología, se puede cambiar el orden en que los grupos son visualizados. 
                  #Por ejemplo, se pueden ordenar las puntas de manera alfabética (en la medida de lo posible), o con los grupos más derivados hacia uno de los lados del árbol. Para escalerizar hacia la derecha:
                  plotTree(ladderize(arbolNJ))
                  
          #GUARDAR UN ARBOL 
            write.tree(arbolNJ, file = "file_name.nex")
          #LEER ARCHIVO
            read.tree(file = "file_name.nex")

#4. METODO CONSTRUCCIÓN -ENRAIZAR
  #Para enraizarlos podemos usar la función root del paquete ape. 
  #Por ejemplo, para poner por raíz las secuencias de fraxatina de Ornitorinco pasamos por argumento el nombre de la secuencia que corresponde a esta secuencia en el parámetro outgroup (Ornitorrinco).
    arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
    plot(arbolNJraiz)
  
  #Se puede hacer lo mismo con el árbol creado a partir del método UPGMA.
  arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
  plot(arbolUPGMAraiz)

    #Visualizar los dos árboles a la vez con los siguientes comandos:
    layout(matrix(c(1,2)), height=c(10,10)) 
    par(mar=c(1,1,1,1))
    plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
    plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)
  
  #4.1 Árboles de Parsimonia
    #Disminuir el número de pasos que explican un árbol evolutivo contando el número de cambios de cada uno de los caracteres. 
    #En los métodos de distancia (UPGMA y NJ) sólo se llega a un único árbol
    # En parsimonia se evalúan múltiples árboles. Se suman todos los pasos reqerridos para explicar un árbol obteniendo y se obtiene un número que se compara entre otros árboles con diferentes pasos. 
    #El árbol con menor número de cambios es el ideal.
      #Hay que tener en cuenta que la parsimonia no usa todos los caracteres (aminoácidos de la secuecia). 
      # Son eliminados aquellos que son constantes en todos los taxones, o aquellos que son variables pero no informativos. 
        #Un caracter es informativo si tiene al menos dos estados de caracter y por lo menos dos de estos estados ocurren con una frecuencia mínima de dos.
  
      #Para estimar árboles de máxima parsimonia existen varias posbilidades, la más sencilla es partir de árboles de distancia. 
      #Se utiliza un árbol de inicio obtenido por distancia y se cuenta su número de pasos. 
      #Podemos estimar, por ejemplo, el número de pasos del árbol arbolUPGMAraiz.
      parsimony(arbolUPGMAraiz, fraxatin) #[1] 313
        #El árbol arbolUPGMAraiz tiene 313 pasos, algo importante es que aunque esté con raíz o no el número de pasos debe ser el mismo. 
        #Probémoslo usando el árbol sin raíz:
        parsimony(arbolUPGMA, fraxatin)  
        
        #ARBOL MEJOR PARSIMONIA 
        
        #1. Método optim.parsimony 
          #Este método permite encontrar árboles bajo máxima parsimonia usando árboles de distancia de inicio.
            mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin)
           #Ahora hagámoslo con el árbol de NJ:
            mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin) # Final p-score 307 after 1 nni operations
           
        #2. Algoritmo de búsqueda pratchet. El cual tiene los siguientes pasos:
            #Generar un árbol de inicio con algún nivel de intercambio de ramas 
            #Seleccionar al azar un subconjunto de caracteres (aminoácidos) y darles más peso, esto quiere decir que se usan con mayor frecuencia dichos caracteres para hacer los análisis. 
            #La cantidad de caracteres que se seleccionan es establecida por el usuario, típicamente entre 5 y 25% de los caracteres. 
            #Realizar intercambio de ramas (branch swapping). Iterar pasos 1 a 3 entre 50 y 200 veces. 
            #El algoritmo anterior es computacionalmente más amigable y rápido. Probémoslo usando la función pratchet:
             fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)
                #Se genera un árbol, se pesan caracteres por defecto y se itera un número determiando de veces. El algoritmo ha encontrado múltiples ocasiones que el árbol más corto es de 307 pasos pero sólo 4 árboles con igual longitud y número de pasos aunque presenten diferente topología.
                fraxatin_parsimonia
                #Para poderlos comparar es necesario enraizarlos.
                fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
                plot(fraxatin_parsimoniaR, cex = 0.6)
                  #Para escoger el sólo un árbol con igual parsimonia hay que conseguir un árbol de consenso. Una posibilidad es un árbol de consenso estricto, sólo los grupos monofiléticos presentes en todos los árboles son incluidos; si un grupo monofilétco no está presente en todos los árboles colapsa formando una politomía en el árbol consenso. Este tipo de consenso es muy riguroso y suele tener muchas politomías, de modo que se puede hacer un consenso con una regla de mayoría. 
                  #Por ejemplo, si un grupo está presente en una mayoría del 51% de árboles, entonces el árbol consenso lo incluye, si sólo está presente en el 50%, el consenso no lo usa.
  
                  #Para hacer un árbol de consenso estricto podemos usar el método ape con parámetro p de 1, que corresponde a un 100% de consenso entre ramas.
                    estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
                    plot(estrictode100, cex = .6)
                    #Para un árbol menos estricto podemos cambiar el valor del parámetro p:
                      estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)
                      plot(estrictode30, cex = .6)
                      #El resultado es un árbol con menos politomías, más resuelto.
  
#5. Bootstrap - Encontrar el mejor árbol y evaluar su fiabilidad. 
    #METODO BASADO EN LA PROBABILIDAD, CREA REPLICAS DE TODOS LOS ARBOLES Y DA UN VALOR 
    #Para dar soporte a los árboles se puede hacer una serie de seudoréplicas con remplazamiento (bootstrapping) de una matriz. Los cambios entre seudoréplicas consisten en el uso diferencial de los caracteres. Se crean árboles en los que un mismo caracter se repite u otro no se usa.
    #Con cada réplica se hace un árbol consenso; la veces que un grupo se repita en el conjunto de réplica es el valor de soporte del nodo. Esta técnica, ampliamente utilizada, proporciona evaluaciones de «confianza» para cada clado de un árbol observado, basándose en la proporción de árboles bootstrap que muestran ese mismo clado. Normalmente un nodo con un soporte mayor a 80% es un buen soporte.
     arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)
      # En este caso usamos la función pratchet, si usáramos otra se demoraría más. También usamos un número de réplicas irrisosio pues este es un ejercicio demostrativo y no es necesario gastar mucho tiempo en él. La rutina anterior genera entonces 10 árboles pseudoréplicas.
       plot(arbolesbootstrap, cex = .6)
      #Ahora bien, generamos un consenso; en este caso con un consenso al 60%:
      estricto60 <- consensus(arbolesbootstrap, p = 0.6)
      plot(estricto60, cex = .6)
  
#6. MODELOS PROBABILISICOS 
    # Métodos de reconstrucción filogenética basados en modelos probabilísticos. 
    #Los caracteres en este caso no son sólo caracteres si no que también representan dinámicas evolutivas. 
    #A diferencia de la parsimonia, con estos modelos todos los caracteres son útiles, incluso los que son constantes entre todos los taxones. Además estos cumplen con un modelo evolutivo.
  
  #6.1 Árboles de máxima verosimilitud
    #Se calcula la verosimilitud (probabilidad de obtener un dato según un modelo) de un árbol de acuerdo a un alineamiento de secuencias usando un modelo de sustitución de aminoácidos. También se tienen en cuenta la frecuencia de de ocurrencia de los aminoácidos. Los pasos para determinar la verosimilitud de un árbol son:
      #Se genera un árbol enraizado de inicio de cualquier tipo, puede ser al azar.
      #Se calcula la verosimilitud de casa sitio (de aminoácidos) usando un árbol.
      #Se calcula la verosimilitud total del árbol por nodo, se consideran todos los escenario posibles que pudieron dar origen a los estados de caracter observados para dicho sitio. La verosimilitud de un sitio es, entonces, la suma de verosimilitudes de la reconstrución de este sitio (en forma de árbol) dado un modelo de evolución o sustitución.
      #La verosimilitud es computacionalmente costosa. Para el ejemplo usamos de nuevo nuestro objeto fraxatin de clase phy y creamos un árbol al azar de 11 ramas (porque tenemos 11 secuencias) con rtree como punto de partida.
      arbolazar <- rtree(n = 11, tip.label = names(fraxatin))
      plot(arbolazar, cex = .5)
        #En seguida lo enraizamos por las secuencias de Ornitorinco para poderlo visualizar mejor. Además los «escalerizamos» hacia la derecha y le agregamos escala; aquí la longitud de la rama sí es significativa, indica cantidad de cambio en cuanto a sustituciones de aminoácidos.
        arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
        plot(ladderize(arbolazarR), cex = .5); add.scale.bar()
        # A partir del arbol de arriba se puede iniciar la búsqueda del mejor árbol por máxima verosimilitud. 
        # Lo primero que se hace es carcular la verosimilitud del árbol dadas las secuencias. Con pml (Phylogenetic maximum likelihood), podemos computar tal verosimilitud.
          ajustado <- pml(arbolazarR, fraxatin)
          ajustado
              #model: Mk 
              ## loglikelihood: -4367.895 
              ## unconstrained loglikelihood: -1479.871 
              ## Rate matrix:
        # La información que tiene el objeto ajustado nos reporta la verosimilitud del árbol al azar que habíamos creado, que es -4348.064. 
        #También reporta un modelo de substitución general, el cual tal vez no se ajuste bien a los datos. 
        #Lo que hay que hacer es encontrar un árbol que optimice la verosimilitud usando un modelo de sustitución; para esto vamos a usar el método optim.pml del paquete phangorn, el cual computa la verosimilitud de un árbol filogenético dado un alineamiento múltiple de secuencias y un modelo de evolución de AA. Toma como argumentos un objeto de clase pml, el tipo de modelo que se quiere usar así como el tiempo de rearreglo para los árboles.
         ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")
            #Para ver el árbol oculto usamos $tree. También lo enraizamos.
            ajustadoconDay$tree
    
           ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")
           plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()
  
  El árbol anterior fue generado usando la matriz de sustitución de Dayhoff. Pero se pueden usar diferentes modelos.
  ajustadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")
  ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")
  Podemos comparar los modelos calculando el Criterio de información de Akaike AIC:
  AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)
  ##                df      AIC
  ## ajustadoconDay 19 5188.886
  ## ajustadoconBlo 19 5197.264
  ## ajustadoconJTT 19 5061.388
  La primera columna corresponde a los grados libertad. Según el criterio anterior, el mejor modelo que se ajusta con los datos (con el AIC más bajo) es JTT modelo de Jones-Taylor-Thornton para evaluar la distancia entre secuencias de proteínas y optimiza la verosimilitud.
  mejorarbol <- optim.pml(
    object = ajustadoconDay, 
    model = "JTT", 
    rearrangement = "ratchet")
  
  mejorarbol
  
  mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
  plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()