# List of the exercises solved in the textbook

The list is written in French, as the textbook. But English translation can be found in the Python documentation.
 
# CHAPTER 1: Equations locales de conservation et d'évolution
 
    1.0 RAPPELS                                                                     {1}
        1.0.1 Equations de conservation                                             {1}
        1.0.2 Equation d'évolution de l'entropie                                    {2}
    1.1 Equation de conservation de l'enthalpie massique h                          {3}
    1.2 Equation d'évolution de l'entropie massique s                               {4}
    1.3 Equations en variables primitives V et T                                    {5}
 
# CHAPTER 2: Modèles mathématiques pour l'aérodynamique
 
    2.0 RAPPELS                                                                     {7}
        2.0.1 Modèles mathématiques                                                 {7}
        2.0.2 Fonction de courant et potentiel des vitesses                         {9}
    2.1 Equations et nombres sans dimension                                         {10}
        2.2  Notion de fluide parfait                                               {12}
        2.3 Fonction de courant et potentiel des vitesses                           {13}
            2.3.1 Tourbillon ponctuel                                               {13}
            2.3.2 Source ponctuelle                                                 {14}
            2.3.3 Doublet source-puits d'axe horizontal                             {15}
            2.3.4 Etirement pur                                                     {16}
        2.4 Atmosphère standard                                                     {17}
 
 # CHAPTER 3: Théorie des potentiels complexes: application aux profils d'aile
 
    3.0 RAPPELS                                                                     {19}
        3.0.1 Potentiel complexe                                                    {19}
        3.0.2 Ecoulements élémentaires                                              {20}
        3.0.3 Profil d'aile                                                         {21}
        3.0.4 Transformation conforme                                               {21}
        3.0.5 Force de portance, moment de tangage                                  {22}
    3.1 Ecoulement autour d'un cylindre et d'une sphère                             {23}
    3.2 Transformation de Joukowski                                                 {26}
    3.3 Transformation de Karman-Trefftz                                            {33}
    3.4 Transformation de Von Mises                                                 {36}
    3.5 Transformation de Van de Vooren et Jong                                     {39}
    3.6 Profil à double pointe et double courbure                                   {42}
    3.7 Plaque plane en incidence                                                   {44}
    3.8 Introduction à la méthode des panneaux                                      {47}
    3.9 Effet de sol                                                                {50}
    3.10 Ailes en tandem                                                            {53}
    3.11 Transformation de Schwarz et Christoffel: expansion dans un canal          {55}
    3.12 Ovale de Rankine                                                           {57}
 
 # CHAPTER 4: Théorie linéarisée des profils minces
 
    4.0 RAPPELS                                                                     {61}
        4.0.1 Introduction                                                          {61}
        4.0.2 Modèle linéarisé                                                      {62}
        4.0.3 Résolution par la méthode de Glauert                                  {63}
        4.0.4 Caractéristiques aérodynamiques                                       {64}
    4.1 Plaque plane en incidence                                                   {65}
    4.2 Profils de Joukowski                                                        {67}
    4.3 Famille de profils minces                                                   {69}
    4.4 Profil mince avec volet et bec                                              {71}
    4.5 Lois de squelette et d'épaisseur                                            {74}
    4.6 Braquage de volet ou de bec sur une aile rectangulaire supposée
        d'envergure infinie                                                         {76}
    4.7 Profil NACA 4 chiffres                                                      {78}
 
 # CHAPTER 5: Théorie de la ligne portante de Lanchester-Prandtl
 
    5.0 RAPPELS                                                                     {81}
        5.0.1 Introduction                                                          {81}
        5.0.2 Théorie de Lanchester-Prandtl                                         {82}
        5.0.3 Méthode de Glauert                                                    {83}
        5.0.4 Loi de Biot-Savart et tube tourbillon}                                {84}
     5.1 Calcul analytique des coefficients A_n                                     {85}
        5.1.1 Aile optimale elliptique                                              {85}
        5.1.2 Vrillage d'une aile elliptique                                        {87}
        5.1.3 Braquage d'un volet sur une aile elliptique                           {92}
        5.1.4 Calcul du foyer sur une aile elliptique                               {95}
     5.2 Sillages                                                                   {97}
        5.2.1 Sillage tourbillonnaire de l'Airbus A380-800                          {97}
        5.2.2 Corrections de soufflerie                                             {99}
        5.2.3 Tourbillon en fer à cheval                                            {102}
        5.2.4 Vol en formation                                                      {106}
        5.2.5 Tourbillons de sillage au voisinage du sol                            {108}
     5.3 Résolution numérique du problème de Prandtl                                {110}
        5.3.1 Cas général d'une aile en milieu infini                               {110}
        5.3.2 Ligne portante : aile avec effet de sol                               {118}

 # CHAPTER 6: Théorie de la surface portante et des corps élancés                  
 
     6.0 RAPPELS                                                                    {121}
         6.0.1 Introduction                                                         {121}
         6.0.2 Théorie de la surface portante                                       {122}
         6.0.3 Théorie des corps élancés                                            {123}
     6.1 Evolution du coefficient de portance d'une aile avec l'allongement         {125}
     6.2 Aile en forme de Delta                                                     {126}
     6.3 Distribution de pression sur un ellipsoïde                                 {127}
 
 # CHAPTER 7: Aspects numériques: méthode des singularités
 
     7.0 RAPPELS                                                                    {133}
         7.0.1 Introduction, méthodologie}{133}
         7.0.2 Potentiel des vitesses}{134}
         7.0.3 Conditions aux limites, approche de Dirichlet ou de Neumann          {134}
     7.1 Profils, solution analytique pour la validation                            {135}
         7.1.1 Profil symétrique de Van de Vooren et Jong                           {135}
         7.1.2 Squelette cambré de forme parabolique                                {138}
     7.2 Etude de différentes distributions de singularités                         {140}
         7.2.1 Méthode des tourbillons discrets: squelette cambré                   {140}
         7.2.2 Distribution de sources ponctuelles,  approche de Neumann: 
               profil de Van de Vooren et Jong sous incidence nulle                 {143}
     7.2.3 Distribution uniforme de sources et de doublets, approche Dirichlet:
               profil de Van de Vooren et Jong en incidence                         {146}
 
 # CHAPTER 8: Ecoulements compressibles subsonique et transsonique
 
     8.0 RAPPELS                                                                    {153}
         8.0.1 Introduction                                                         {153}
         8.0.2 Régime compressible subsonique 0<M_0<0.8                             {154}
         8.0.3 Régime transsonique M_0>0.8                                          {154}
         8.0.4 Aile d'envergure finie et aile en flèche                             {154}
     8.1 Ecoulement subsonique isentropique                                         {155}
         8.1.1 Pression totale et pression isentropique                             {155}
         8.1.2 Tube de Pitot en régime compressible                                 {156}
     8.2 Ecoulement transsonique                                                    {158}
         8.2.1 Corrections de compressibilité                                       {158}
         8.2.2 Mach critique inférieur : effet de flèche                            {160}
     8.3 Ecoulement au voisinage d'une paroi ondulée en régime subsonique           {162}
     8.4 Ecoulement au voisinage d'une paroi ondulée en régime transsonique         {164}
     8.5 Ecoulement transsonique autour d'un profil                                 {171}
 
 # CHAPITRE 9: Ecoulements supersoniques linéarisés
 
     9.0 RAPPELS                                                                    {173}
         9.0.1 Introduction                                                         {173}
         9.0.2 Relations fondamentales et méthodes des caractéristiques             {173}
         9.0.3 Applications aux profils                                             {175}
     9.1 Applications aux écoulements externes                                      {176}
         9.1.1 Plaque plane                                                         {176}
         9.1.2 Ecoulement proche d'une paroi ondulée : régime supersonique          {177}
         9.1.3 Caractéristiques d'un profil losangique                              {179}
         9.1.4 Profil de traînée minimale                                           {182}
         9.1.5 Profil lenticulaire                                                  {186}
         9.1.6 Interaction jet supersonique-profil                                  {187}
     9.2 Ecoulements en canal plan                                                  {191}
         9.2.1 Plaque plane en incidence dans un canal supersonique                 {191}
         9.2.2 Canal supersonique divergent                                         {193}
         9.2.3 Profil losangique placé dans un canal supersonique                   {194}
         9.2.4 Remarque sur le plan de l'hodographe                                 {196}
 
 # CHAPTER 10: Ecoulements compressibles monodimensionnels
 
     10.0 RAPPELS                                                                   {197}
         10.0.1 Introduction                                                        {197}
         10.0.2 Lois isentropiques                                                  {198}
         10.0.3 Relations du choc droit                                             {198}
         10.0.4 Ecoulement monodimensionnel dans une tuyère                         {198}
         10.0.5 Ecoulement monodimensionnel avec frottement : lois de Fanno         {199}
         10.0.6 Ecoulement monodimensionnel avec apport de chaleur: 
                lois de Rayleigh                                                    {201}
     10.1 Calcul d'un choc droit                                                    {202}
         10.1.1 Exemple élémentaire                                                 {202}
         10.1.2 Compression par choc et compression isentropique                    {202}
    10.2 Ecoulement dans une conduite                                               {203}
         10.2.1 Ecoulement de Fanno : régime amont subsonique                       {203}
         10.2.2 ecoulement de Fanno : régime amont supersonique                     {205}
         10.2.3 Ecoulement de Rayleigh                                              {206}
         10.2.4 Ecoulement de Rayleigh : apport de chaleur par combustion           {207}
     10.3 Mesure de vitesse par tube de Pitot                                       {209}
         10.3.1 Régime subsonique                                                   {209}
         10.3.2 Régime supersonique                                                 {210}
     10.4 Tuyères                                                                   {211}
         10.4.1 Différents régimes d'écoulement dans une tuyère                     {211}
         10.4.2 Poussée d'une tuyère de propulseur                                  {214}
         10.4.3 Tuyère supersonique et tube de Pitot                                {219}
         10.4.4 Moteur de fusée                                                     {220}
 
 # CHAPTER 11 : Ecoulements supersoniques bidimensionnels
 
     11.0 RAPPELS                                                                   {221}
         11.0.1 Introduction                                                        {221}
         11.0.2 Choc oblique                                                        {222}
         11.0.3 Détente de Prandtl-Meyer                                            {223}
         11.0.4 Coefficients aérodynamiques                                         {223}
         11.0.5 Diagramme pression-déviation                                        {224}
     11.1 Ecoulement autour d'un profil                                             {225}
         11.1.1 Coefficients aérodynamiques d'un profil en incidence                {225}
         11.1.2 Plaque plane en incidence                                           {226}
         11.1.3 Ecoulement supersonique autour d'un profil losangique               {229}
         11.1.4 Méthode choc-détente pour un profil lenticulaire                    {233}
     11.2 Autres géométries                                                         {235}
         11.2.1 Ecoulement supersonique dans un canal: réflexions d'onde de choc    {235}
         11.2.2 Interaction de deux chocs obliques de la même famille               {237}
         11.2.3 Marche descendante en écoulement supersonique                       {239}
     11.3 Réflexion et interactions                                                 {241}
         11.3.1 Réflexion d'une onde de choc sur une ligne isobare                  {241}
         11.3.2 Réflexion d'un faisceau de détente sur une ligne isobare            {242}
         11.3.3 Réflexion d'un faisceau de détente sur une paroi solide             {244}
         11.3.4 Interaction onde de choc oblique - ligne de glissement              {246}
     11.4 Tuyères                                                                   {249}
         11.4.1 Tuyère : interaction d'ondes de chocs en sortie                     {249}
         11.4.2 Tuyère : interaction d'un faisceau de détente en sortie             {251}
         11.4.3 Tracé de lignes caractéristiques dans le plan de l'hodographe       {253}
 
 # CHAPTER 12: Méthodes des caractéristiques en régimes stationnaire et instationnaire
 
     12.0 RAPPELS                                                                   {255}
          12.0.1 Introduction                                                       {255}
          12.0.2 Ecoulement stationnaire bidimensionnel                             {256}
          12.0.3 Ecoulement instationnaire monodimensionnel                         {257}
     12.1 Ecoulement supersonique dans un canal plan                                {257}
     12.2 Tuyère de longueur minimale                                               {262}
     12.3 Piston uniformément accéléré                                              {268}
     12.4 Tube à choc                                                               {271}
     
 # CHAPTER 13: Corps élancés en écoulements supersoniques
 
     13.0 RAPPELS                                                                   {275}
         13.0.1 Introduction                                                        {275}
         13.0.2 Corps élancés: potentiel et coefficient de pression                 {276}
         13.0.3 Optimisation de la traînée sur un corps élancé                      {277}
     13.1 Solutions élémentaires pour un point source                               {278}
     13.2 Obstacle conique                                                          {280}
     13.3 Ogive parabolique                                                         {283}
     13.4 Optimisation sur le corps de Sears                                        {286}
     13.5 Mach amont minimal sur un corps conique                                   {291}
 
 # CHAPTER 14: Ecoulements hypersoniques
 
     14.0 RAPPELS                                                                   {293}
         14.0.1 Introduction                                                        {293}
         14.0.2 Onde de choc droite                                                 {293}
         14.0.3 Onde de choc oblique et similitude hypersonique                     {294}
         14.0.4 Détente de Prandtl-Meyer                                            {295}
         14.0.5 Choc conique en hypersonique                                        {296}
         14.0.6 Méthode de Newton                                                   {297}
     14.0.7 Méthodes de calcul                                                      {297}
     14.1 Profils losangique et triangulaire                                        {298}
     14.2 Plaque plane                                                              {301}
     14.3 Ecoulement autour de cones ou de dièdres                                  {303}
     14.4 Méthodes de calcul pour les obstacles pointus                             {305}
     14.5 Traînée hypersonique minimale                                             {309}
     14.6 Méthodes de calcul pour les obstacles émoussés                            {312}
 
 # CHAPTER 15: Effets visqueux et couche limite
 
     15.0 RAPPELS                                                                   {317}
         15.0.1 Introduction                                                        {317}
         15.0.2 Régime laminaire                                                    {318}
         15.0.3 Régime transitionnel                                                {319}
               1. Critère de Michel                                                 {320}
               2. Critère de Cebeci et Smith                                        {320}
               3. Critère de Granville                                              {321}
         15.0.4 Régime turbulent                                                    {321}
     15.1 Régime laminaire                                                          {324}
         15.1.1 Plaque plane: loi approchée en sinus                                {324}
         15.1.2 Plaque plane : loi polynomiale                                      {325}
         15.1.3 Aspiration uniforme d'une couche limite de plaque plane             {329}
         15.1.4 Mise en mouvement brusque d'une plaque plane                        {330}
         15.1.5 Plaque plane : relations intégrales avec soufflage/aspiration       {332}
         15.1.6 Ecoulement de Couette instationnaire                                {334}
         15.1.7 Ecoulement oscillant : deuxième problème de Stokes                  {337}
         15.1.8 Solution de similitude avec gradient de pression : 
                solution de Falkner et Skan                                         {339}
         17.1.9 Ecoulement laminaire dans un jet                                    {343}
     15.2 Abscisse de transition dans une couche limite de plaque plane             {347}
     15.3 Régime turbulent                                                          {349}
         15.3.1 Profil de vitesse turbulent: loi en puissance                       {349}
         15.3.2 Profil de vitesse turbulent, loi de Van Driest                      {353}
         15.3.3 Fluctuations turbulente en proche paroi                             {355}
         15.3.4 Turbulence : soufflage avec gradient de pression                    {357}
         15.3.5 Couche limite turbulente avec aspiration                            {358}
         15.3.6 Ecoulement turbulent dans un jet                                    {362}
 
  
 
