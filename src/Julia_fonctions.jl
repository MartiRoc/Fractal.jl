using Images
using Dates

######################### FONCTIONS GENERATRICES DE FRACTALES #########################

######################### par ITERATIONS INVERSES, c.f. la partie théorique du rapport

#=
Cette première fonction sert à calculer le point de départ de l'algorithme d'itérations
inverses qui approche un ensemble de Julia de paramètre c, il s'agit de partir du point 
fixe du polynôme -- X -> X^2 + c -- qui admet la plus grande valeur du module de la dérivée 
-- X -> 2X --. C'est judicieux pour diverses raisons théoriques, la première étant que ce 
point appartient forcément à la fractale, ce qui est une condition nécéssaire du 
théorème sur lequel se base l'aglorithme. 

Paramètre(s):
- c : complexe, paramètre de la fractale.
=#

function PointDeDepart(c)
    if c == 0.25
        return 0.5
    else
        d = sqrt(1-4*c)
        z1 = (1-d)/2
        z2 = (1+d)/2
        if abs(2*z1) > abs(2*z2)
            return z1
        else
            return z2
        end
    end
end

#= 
1-er G : cette fonction descend l'arbre binaire des pré-images à la manière du pachinko, 
on choisit aléatoirement une branche, N désigne donc la profondeur à laquelle on va 
descendre mais aussi le nombre de points calculés. Cette fonction approche l'ensemble de 
julia de paramètre c.

Paramètre(s):
- c : complexe, paramètre de la fractale,
- N : entier, nombre de points calculés,
- z0 : complexe, point de départ
=#

function Julia1(N, c ; z0 = -42)

    # A moins que z0 ne soit renseigné, il est calculé :
    (z0 == -42) && (z0 = PointDeDepart(c)) 

    # Déclaration des vecteurs X -> partie réelle, Y -> partie imaginaire des points calculés.
    X = [real(z0);zeros(N-1)]
    Y = [imag(z0);zeros(N-1)]

    # Descente de l'arbre des pré-images.
    for i in 2:N
        z0 = ((-1)^rand(0:1))*sqrt(z0-c)
        X[i] = real(z0)
        Y[i] = imag(z0)
        #(mod(i,N/10)==0) && print(string((i/N)*100),"%","\n")
    end

    return [X,Y]
end

#=
2-ième G : Après quelques générations d'images, descendre l'arbre des inverses uniquement 
en profondeur (1-er G) pose un souci, certaines "zones" de la fractales semblent ne pas 
apparaître, il faut aussi parcourir l'arbre en largeur ! La fonction suivante calcule TOUTES 
les pré-images jusqu'à une profondeur P. Attention on calcul donc 2^P+1 - 1 points. 

Paramètre(s):
- c : complexe, paramètre de la fractale,
- P : entier, profondeur de la descente, 2^P+1 - 1 points calculés.
=#

function Julia2(P, c) 

    z0 = PointDeDepart(c)
    X = [[real(z0)];zeros(2^(P+1)-2)]
    Y = [[imag(z0)];zeros(2^(P+1)-2)]

    for d in 1:P
        for i in 2^(d-1):2^d-1
            z2 = sqrt(X[i]+Y[i]im-c)
            X[2*i] = real(z2)
            Y[2*i] = imag(z2)
            X[2*i+1] = real(-z2)
            Y[2*i+1] = imag(-z2)
            
        end
        #print(string((d/P)*100),"%","\n")
    end

    return [X,Y]
end

#=
3-ième G : Le précédent algorithme descend l'arbre des inverses intégralement jusqu'à 
une profondeur P, mais la progression exponentielle du nombre de points calculés nous
empêche de descendre "profondemment" dans l'arbre. Ce troisième G fait un compromis
entre une descente en profondeur (Julia1) et en largeur (Julia2). L'idée est de 
descendre à la manière de Julia1 à une profondeur N, mais à chaque point calculé on
relance Julia1 jusqu'à une profondeur N. On réitère l'opération B fois. On descends ainsi 
jusqu'à une profondeur 2N, et le caractère aléatoire de Julia1 qui est lancé B*N*N fois nous 
assure de couvrir "une bonne largeur". Le nombre de points différents calculés est 
majoré par B*N^2 (majoré car les "descentes" vont se croiser, surtout au début de l'arbre).

Paramètre(s):
- c : complexe, paramètre de la fractale,
- N , B : entiers, paramètres de la descente de l'arbre des pré-images, B*N^2 points calculés.
=#

function Julia3(N, c; B = 5)

    z0 = PointDeDepart(c)
    X = zeros(B*N*N)
    Y = zeros(B*N*N)

    for b in 1:B
        Temp1 = Julia1(N,c)
        rang = (b-1)*N*N
        for i in 1:N
            Temp2 = Julia1(N,c, z0 = Temp1[1][i]+Temp1[2][i]im)
            X[rang + (i-1)*N + 1 : rang + (i-1)*N + N] = Temp2[1]
            Y[rang + (i-1)*N + 1 : rang + (i-1)*N + N] = Temp2[2]
        end
    end

    return [X,Y]
end

######################### TESTS (environ même nombre de points calculés : 2.1M)

#M = Gray.(VecToMat(Julia1(2100000,-0.66+0.4im)))
#M = Gray.(VecToMat(Julia2(20,-0.66+0.4im)))
#M = Gray.(VecToMat(Julia3(648,-0.66+0.4im)))

#=
Remarque(s): C'est rapide, pas compliqué à implémenter, mais pas fameux quant aux détails qu'on
arrive à "capter" de la fractale. Même en poussant le nombre de points calculés (en profondeur
et en largeur), les points "centraux" de la fractale semblent inaccessible en un temps raisonnable.
L'algorithme suivant donne de biens meilleurs résultats. Spoiler : 
=#

#M = Gray.(VecToMat(Julia4Rv2(-0.66+0.4im)))

######################### par TEMPS DE FUITE, c.f. la partie théorique du rapport

#=
Cet algorithme tire partie directement de la définition de l'ensemble de julia de paramètre c.
En effet c'est la frontière de l'ensemble des complexes pour lesquels les itérations du poly-
nôme -- X -> X^2 + c -- restes bornées. L'idée de l'algorithme et de quadriller le plan complexe
et de tester pour chacun des éléments du quadrillage si itérérer I fois (à choisir) le polynôme
ci-dessus conduit à les faire sortir d'un cercle de rayon bien choisi (2), impliquant la divergence
si nous continuyons à itérer. On approche ainsi l'ensemble de julia REMPLI de paramètre c par
"élimination". Les éléments qui ont divergé ne sont pas dans l'ensemble de julia, ceux qui restent
sont soit dans l'ensemble de julia rempli, soit ne le sont pas mais ont "réussi" le test des I 
itérations grâce à leur vitesse (lente) de divergence. Le choix de I est donc directement lié à la 
précision de l'approximation que nous obtenons.
Pour passer ensuite de l'ensemble de julia rempli approché à l'ensemble de julia (approché), il 
suffit de prendre la frontière de l'ensemble des points retenus. La fonction suivante fournit des
points appartenant à l'ensemble de julia rempli.

Paramètre(s):
- c : complexe, paramètre de la fractale,
- pas : réel, pas du quadrillage de [-2,2]x[0,2],
- I : entier, nombre d'itérations maximal pour qu'un point du quadrillage réussisse "le test de 
divergence".
=#

function Julia4R(c ; pas = 0.0004, I = 50)

    # Le quadrillage :
    Px = collect(range(-2,2,step=pas))
    Py = collect(range(0,2,step=pas))
    #= Remarque : on ne quadrille que le pan supérieur, on obtient le reste par symétrie 
    (si z est dans l'ensemble, -z l'est aussi par définition 
    et puisque z*z + c = (-z)*(-z) + c). =#

    #= Initialisation des vecteurs résultats, X -> partie réelle, Y -> partie imaginaire 
    des points calculés : =#
    N = size(Px)[1]*size(Py)[1]
    X = zeros(N).+0.0
    Y = zeros(N).+0.0

    r = 1 # Un compteur qui garde en mémoire jusqu'oû seront remplis les vecteurs ci-dessus.
    cx = real(c)
    cy = imag(c)

    for x in Px
        for y in Py
           iteration = 0
           x0 = x
           y0 = y
           x1 = x*x
           y1 = y*y
           while ((x1 + y1 <= 4) & (iteration < I))
                #= On utilise ici des petits résultats sur les complexes pour minimiser le 
                nombre d'opérations, surtout de multiplications. =#
                y0 = (x0+x0)*y0+cy
                x0 = x1-y1 + cx
                x1 = x0*x0
                y1 = y0*y0
                iteration += 1
           end
           if iteration == I # Le point ne s'est pas "échappé" après I itérations.
                X[r]=x
                Y[r]=y
                r+=1
           end
        end    
    end

    # On applique la symétrie.
    X = [X[1:r-1] ; -X[1:r-1]]
    Y = [Y[1:r-1] ; -Y[1:r-1]]

    return [X ,Y]
end

#=
Pour être plus précis, mieux vaut utiliser la fonction ci-dessous que diminuer le pas de 
la fonction ci-dessus (ce qui conduirait à l'évaluation de beaucoup de points inutiles). 
La fonction ci-dessous restreint la région du plan evaluée à un rectangle proche de la 
fractale en appelant Julia1 au préalable. On a un coup en ressource "fixe" plus important 
mais on est gagnant sur un quadrillage plus fin du plan.

Par défaut Julia4R quadrille [-2,2]x[0,2] en 10000x5000 = 50 000 000 de points, tandis que par 
défaut Julia4Rv2 quadrille un rectangle strictement inclu dans [-2,2]x[0,2] en 10000x5000 de 
points. Julia4Rv2 est donc par défaut "plus fine" que Julia4R. L'intêrét de Julia4Rv2 se 
révèle lorsque par exemple je divise le pas dans Julia4R par 10 (=0.00004), l'algorithme 
n'aboutit pas à cause d'un problème de mémoire (du à ma configuration) alors qu'augmenter la 
finesse du quadrillage d'un facteur 5 dans Julia4Rv2 suffit à obtenir le même niveau de détail 
que si Julia4R avait abouti (car la plupart des fractales ne s'étendent que sur la moitié en 
hauteur ou largeur du rectangle [-2,2]x[-2,2]) et lui aboutit.

Paramètre(s):
- c : complexe, paramètre de la fractale,
- len : entier, dimension du quadrillage = len*len*0.5,
- I : entier, nombre d'itérations maximal pour qu'un point du quadrillage réussisse "le test de 
divergence".
=#

function Julia4Rv2(c ; len = 10000, I = 50)

    #= On restreint le quadrillage à un quadrillage proche de la fractale en appelant Julia1
    qui est peu gourmande en ressource. =#
    Approx = Julia1(1000, c)
    minX = minimum(Approx[1])-0.05
    maxX = maximum(Approx[1])+0.05
    maxY = maximum(Approx[2])+0.05

    # Quadrillage d'un rectangle proche de la fractale (pan supérieur).
    Px = collect(range(minX,maxX,length=len))
    Py = collect(range(0,maxY,length=Int(round(len*0.5))))

    # Nombre de points dont on va tester la divergence. 
    N = len * Int(round(len*0.5))

    X = zeros(N).+42.0
    Y = zeros(N).+42.0

    r = 1 # Un compteur qui garde en mémoire jusqu'oû seront remplis les vecteurs ci-dessus.
    cx = real(c)
    cy = imag(c)

    for x in Px
        for y in Py
           iteration = 0
           x0 = x
           y0 = y
           x1 = x*x
           y1 = y*y
           while ((x1 + y1 <= 4) & (iteration < I))
                y0 = (x0+x0)*y0+cy
                x0 = x1-y1 + cx
                x1 = x0*x0
                y1 = y0*y0
                iteration += 1
           end
           if iteration == I
                X[r]=x
                Y[r]=y
                r+=1
           end
        end    
    end

    X = [X[1:r-1] ; -X[1:r-1]]
    Y = [Y[1:r-1] ; -Y[1:r-1]]

    return [X,Y]
end

# (***) c = 0.59 + 0.4im
#M = Gray.(VecToMat(Julia4Rv2(0.59+0.4im, I=50)))
#M = Gray.(VecToMat(Julia4Rv2(0.59+0.4im, I=20)))
#M = Gray.(VecToMat(Julia3(648, 0.59+0.4im)))
# (***)

#M = Gray.(VecToMat(Julia4Rv2(0.99+0im, I=10)))
#M = Gray.(VecToMat(Julia3(648, 0.99+0im)))

#M = Gray.(VecToMat(Julia4Rv2(0.39+1im, I=15)))
#M = Gray.(VecToMat(Julia3(648,0.39+1im)))

#M = Gray.(VecToMat(Julia4Rv2(0+0.8im, I=50)))
#M = Gray.(VecToMat(Julia3(648, 0+0.8im)))

#= 
PROBLEME de Julia4R et Julia4Rv2 : si la fractale "recouvre" peu d'espace (***), les points 
s'échappent très vite autour des points "peu nombreux" de la fractale et il est difficile d'en 
faire une représentation approchée par cet algorithme sans "faciliter" le test de divergence 
i.e abaisser I ou beaucoup augmenter la finesse du quadrillage (ce qui moralement conduit à 
tester "plus de points plus proches de la fractales" qui divergent plus lentement). 
Par exemple pour (***), I=50 conduit à la génération d'une image noire. 
=#

######################### GENERATION D'IMAGES #########################
 
# Chemin vers lequel stocker les images

#PATH = "C:\\Cours\\UGA22-23\\LR_projet S1\\Images\\"

######################### MATRICE DE PIXELS

#=
REMARQUE:
Nous voulons à partir des vecteurs X et Y générés par les fonctions Julia1, 2, 3 et 4
générer une image. Il s'agit de tracer dans le plan les points de coordonnées (X,Y). 
Toutefois pour que les fractales soient riches en détails nous calculons un grand nombre 
de pré-images/ points (beaucoup plus que le nombre de pixel / la définition des écrans 
standards), beaucoup de ces pré-images/ points vont correspondre au même pixel à l'écran 
(princippe des tiroirs). Après quelques tentatives infructueuses (chronophages) il semble 
donc judicieux de ne pas tracer tous les points (X,Y) tel que le ferait la fonction plot 
ou scatter du package Plots. 

L'idée est de créer une matrice de pixels à partir des vecteurs X et Y PUIS d'interpréter 
cette matrice comme une image. La fonction suivante construit à partir des vecteurs de 
coordonnées X et Y générés par un des algorithmes ci-dessus la matrice de pixels voulue.

Paramètre(s):
- V : vecteur de deux vecteurs X et Y de coordonnées réelles,
- xa et ya : délimitent le rectangle du plan que représente la matrice, resp. les abscisses 
et les ordonnées (permet de "rogner"),
- xc et yc : sont des coefficients de dilatation de la fractale, xc > 1 étirera la fractale 
en abscisse,
- fl et fL : sont les idicateurs de format, par défaut égaux à 10, mettre fl = 16 et fL = 9 
fournira une image en 16:9ième. ATTENTION, cela déforme la fractale (le repère n'est plus normé),
- r : est le coefficient de définition de l'image, le nombre de pixels est 100*fl*r horizontal, 
100*fL*r vertical. 
=#

function VecToMat(V ; xa=(-2,2), ya=(-2,2), xc=1, yc=1, fl=10, fL=10, r=2)

    # Dilation de la fractale.
    X = V[1]*xc
    Y = V[2]*yc

    #= On discrimine les points de X et Y qui sortent de xa*ya pour éviter des soucis
    d'indices plus bas lors de la conversion des points (x,y) en pixels de la matrice. =#
    b = size(V[1])[1]
    for i in 1:b
        x = X[i]
        y = Y[i]
        if (x>xa[2])|(x<xa[1])|(y>ya[2])|(y<ya[1])
            X[i]=42.
        end
    end

    # Définition de la matrice de pixels qui représente xa*ya :
    rapport_xa_ya = (xa[2]-xa[1])/(ya[2]-ya[1]) #= Pour que la matrice représente bien 
    le rectangle xa*ya sans déformation. =# 
    l = Int(round(100*fl*r*rapport_xa_ya)) # Largeur / nombre de pixels verticaux.
    L = Int(round(100*fL*r)) # Longueur / ... horizontaux.
    M = zeros(L,l)*1.0

    #= Parcours simultané des vecteurs X et Y. On convertit les points (x,y) en coordonnées 
    entières à l'échelle de la matrice. Comme on fait des arrondis, pour éviter les pixels (0,.)
    ou (.,0) qui n'existent pas dans la matrice (indices qui commencent à 1) on utilise max. =#
    for i in 1:b
        if X[i]!=42
            a = L - Int(round((Y[i]-ya[1])*L/(ya[2]-ya[1]))) 
            b = Int(round((X[i]-xa[1])*l/(xa[2]-xa[1])))
            M[max(1,a),max(1,b)] = 1. 
        end
    end

    return M
end

######################### TESTS 

#using BenchmarkTools
#@benchmark M = Gray.(VecToMat(Julia1(10000,0.6+0.4im)))
#@benchmark img = Gray.(VecToMat(Julia4R(-0.7+0.388im), r = 6))
#@benchmark img = Gray.(VecToMat(Julia4Rv2(-0.7+0.388im), r = 6))
#save(PATH*"julia_"*string(Dates.day(now()),"-",Dates.month(now()),"-",Dates.hour(now()),"h",Dates.minute(now()))*".png",img)

######################### GENERATION / COLORISATION

#=
L'idée est simple, on prend une matrice de pixels 0/1 réels, on la transforme en quelque chose que 
connait le package Images suivant la colorisation choisie puis on renvoie l'image. 

Paramètre(s):
- M : matrice de 0/1 à représenter,
- bg : HSV(Hue (réel[360]), Saturation (0-1), value (0-1)), triplet qui correspond à la couleur
du fond (les 0 de la matrice), par défaut noir, 
- fg : HSV(.,.,.), couleur de "l'image" (les 1 de la matrice), par défaut blanc,
- rainbow : booléen, irisation horizontale (arc-en-ciel) de l'image (les 1 de la matrice), 
selon les paramètres h, s, v et a. Par défaut les couleurs iront du rouge (extrême gauche) au 
rouge (extrême droite) en passant par toutes les couleurs une fois,
- h : réel[360], couleur de départ de l'arc en ciel (extrême gauche),
- s : réel 0-1, saturation des couleurs de l'arc-en-ciel, 
- v : réel 0-1, valeur des couleurs de l'arc-en-ciel,
- a : réel, déplacement angulaire sur le cercle HSV, l'arc-en-ciel prendra toute les couleurs
de l'angle h à l'angle h+a (possibilité de faire plus d'un tour).
=#

function MatToImage(M ; bg = HSV(0,0,0), fg = HSV(0,0,1), rainbow = false, h=0, s=1, v=1, a=360)

    # Création de l'image finale.
    img = RGB.(M)

    if !rainbow # Pas d'arc-en-ciel :'(.

        for i in axes(M, 1)
            for j in axes(M, 2)
                # Si le pixel est blanc (1), colorier en "fg".
                if M[i, j] == 1
                    img[i, j] = fg
                # Si le pixel est noir (0), colorier en "bg".
                else
                    img[i, j] = bg
                end
            end
        end

        return img

    else # ARC-EN-CIEL

        #= On cherche le nombre de colonnes de M possèdant un pixel 1 (pour créer la palette
        de couleurs irisée. =#
        l = 0
        for c in axes(M,2)
            if 1 in M[:,c]
                l += 1
            end
        end

        # Si il n'y a pas de 1 dans M, on renvoie l'image noire (que des 0).
        (l == 0)&&(return img)

        # La palette contenant l'arc-en-ciel. 
        palette = range(HSV(h,s,v), stop=HSV(h+a,s,v), length=l)
        
        # Va nous servir à changer de couleur après chaque colonne possédant un 1. 
        compt = 1 

        for i in axes(M, 2) # Parcours des colonnes. 
            bool = false # A-t-on colorié un pixel dans cette colonne ? 
            for j in axes(M, 1) # Parcours des lignes. 
                if M[j, i] == 1
                    img[j, i] = palette[compt]
                    bool = true # On a colorié un pixel dans cette colonne. 
                else
                    img[j, i] = bg
                end
            end
            # Si on a colorié un pixel dans cette colonne, on changera de couleur dans la palette.
            if bool
                compt += 1 
            end     
        end

        return img
    end
end

######################### TESTS 

#img = MatToImage(VecToMat(Julia4Rv2(-0.66+0.4im, len = 10000, I = 100), r = 4, fL = 9, fl = 16, yc = 1.5),bg = HSV(0,0,0.1), rainbow = true, h = 60, a=360, s = 1)
#img = MatToImage(VecToMat(Julia4Rv2(-0.66+0.4im, len = 10000), r = 4),bg = HSV(0,0,0.1), rainbow=true)
#save(PATH*"julia_"*string(Dates.day(now()),"-",Dates.month(now()),"-",Dates.hour(now()),"h",Dates.minute(now()))*".png",img)

######################### FONCTION FINALE

#=
Paramètres :
- c : complexe, paramètre de la fractale voulue,
- algo : "IV" ou "ET", l'algorithme voulu, resp Julia3 (Itérations Inverses) ou 
Julia4Rv2 (Escape Time), par défaut le dernier,
- tous les paramètres par défaut de Julia3 (+N), Julia4Rv2, VecToMat, MatToImage,
- dl : booleen -> sauvegarder l'image sur le PC ? Par défaut "false",
- PATH : chemin vers lequel stoquer l'image, par défaut "", l'image est sauvegardée dans 
le répertoire de travail ou dans C:\Users\*vous* sur windows si aucun répertoire de travail n'a
été spécifié,
- tips : booléen -> voulez-vous des conseils sur l'algorithme "ET" ?
=#

function FractaleMR(c ; algo = "ET", dl = false, PATH = "",
    N = 648, B = 5, #Julia3
    len = 10000, I = 25, #Julia4Rv2
    xa=(-2,2), ya=(-2,2), xc=1, yc=1, fl=10, fL=10, r=2, #VecToMat
    bg = HSV(0,0,0), fg = HSV(0,0,1), rainbow = false, h=0, s=1, v=1, a=360, #MatToImage
    tips = true) 

    if tips
        print("\n", 
        
        "CONSEILS POUR L'ALGORITHME ESCAPE TIME:", "\n", "\n",

        "- Si vous choisissez algo = 'ET' et que la fractale n'apparait", "\n",
        "  pas, tenter de baisser I < 50 (ou augmenter len > 10000).", "\n","\n",

        "- Si vous voyez des zones plaines, essayez d'augmenter I > 50 (ou len).", "\n", "\n",

        "- Si des lignes de la couleur du fond aparaissent dans les zones", "\n",
        "  plaines, c'est que, pour une même surface, la finesse du", "\n",
        "  quadrillage pour le calcul des points -paramètre len- est", "\n",
        "  inférieure à la définition de l'image -paramètre r-, augmentez le", "\n",
        "  premier ou baisser le second.", "\n", "\n",
        
        "- Utiliser -rainbow = true- !", "\n", "\n")
    end

    if algo == "IV"
        img = MatToImage(VecToMat(Julia3(N,c,B=B), 
        xa=xa, ya=ya, xc=xc, yc=yc, fl=fl, fL=fL, r=r),
        bg = bg, fg = fg, rainbow = rainbow, h=h, s=s, v=v, a=a)
    else
        img = MatToImage(VecToMat(Julia4Rv2(c, len=len, I=I), 
        xa=xa, ya=ya, xc=xc, yc=yc, fl=fl, fL=fL, r=r),
        bg = bg, fg = fg, rainbow = rainbow, h=h, s=s, v=v, a=a)
    end

    if dl == true
        save(PATH*"julia_"*string(Dates.day(now()),"-",Dates.month(now()),"-",Dates.hour(now()),"h",Dates.minute(now()))*"_"*"c="*string(real(c))*"+i"*string(imag(c))*".png",img)
    end

    return img
end

######################### TESTS 

# paramètres de Julia3

#FractaleMR(0.99+1im, algo = "IV")
#FractaleMR(0.99+1im, B = 10, N = 648, algo = "IV")

# paramètres de Julia4Rv2

#FractaleMR(0.19+0.6im)
#FractaleMR(0.19+0.6im, I = 100)
#FractaleMR(0.19+0.6im, I = 25)
#FractaleMR(0.19+0.6im, len = 10, I = 100)
#FractaleMR(0.19+0.6im, len = 5000, I = 100)
#FractaleMR(0.19+0.6im, len = 2000, I = 100)

# paramètres de VecToMat

#FractaleMR(0.19+0.6im, I = 100)
#FractaleMR(0.19+0.6im, I = 100, yc = 2)
#FractaleMR(0.19+0.6im, I = 100, xc = 0.1)
#FractaleMR(0.19+0.6im, I = 100, xa = (-1,1), ya=(-1,1))
#FractaleMR(0.19+0.6im, I = 100, fl = 16, fL=9)
#FractaleMR(0.19+0.6im, I = 100, r = 0.1)

# paramètres de MatToImage
#FractaleMR(0.19+0.6im, I = 100, bg=HSV(190,1,0.3), fg=HSV(60,1,1), tips = false)
#FractaleMR(0.19+0.6im, I = 100, xc = 1.5, yc=1.5, rainbow=true, tips=false)
#FractaleMR(0.19+0.6im, I = 100, xa = (-1,1), rainbow=true, h=120, a=120, tips=false)

# sauvegarde
#FractaleMR(0.19+0.6im, I = 100, dl=true, PATH = "C:\\Cours\\UGA22-23\\LR_projet S1\\Images\\")