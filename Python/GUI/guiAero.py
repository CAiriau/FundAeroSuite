
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 
import os
from sys import *
import platform
import numpy             as np
import matplotlib.pyplot as plt


#import tkinter
from tkinter.messagebox import * 
from tkinter.ttk import * 
from tkinter import filedialog  
from tkinter import * 
from contextlib import contextmanager


from Tools.misc                 import *
from CompressibleFlow.fonctions import *
from Livre.Chapter2             import *
from Livre.Chapter3             import *
from Livre.Chapter4             import *
from Livre.Chapter5             import *
from Livre.Chapter6             import *
from Livre.Chapter7             import *
from Livre.Chapter7c            import *
from Livre.Chapter8             import *
from Livre.Chapter9             import *
from Livre.Chapter10            import *
from Livre.Chapter11            import *
from Livre.Chapter12            import *
from Livre.Chapter13            import *
from Livre.Chapter14            import *
from Livre.Corrections          import *
from Livre.misc                 import *
from Choc_Conique.choc_conique  import *
from Livre.panneaux_intro       import *
from Livre.Fortran.fortran      import *
import CheminsParDefaut



#*******************************************
def set_params_gui():
    #*******************************************
    """
    Set default parameters

    """
    prm={}  # Dictionnary of the parameters  

    # physical parameters

    prm["SrcDir"]         = FortranPath                 # chemin pour les sources fortran
    prm["RunDir"]         = FortranRun                  # chemin pour lancer la solution en fortran
    prm["PyRun"]          = PyRun
    
    prm["ShowMenu"]       = True                        # Utilisation d'un menu
    prm["ShowTree"]       = True                        # Utilisation d'un tableau
    prm["ShowLang"]       = False                       # inutile
    prm["ShowFort"]       = True                        # Utilisation du fortran
    prm["ShowCheck"]      = False                       # Utilisation de case  à cocher
    prm["ShowText"]       = False                       # Position d'un texte (aide) dans la fenêtre principale
    prm["ShowTitre"]      = False                       # Affichage du titre du livre
    prm["ShowDefaultCase"] = True                       # Affichage du cas par défaut
    prm["geometry"]       = '1000x800+100+50'           # taille de la fenêtre principale
    prm["redirection"]    = False
    prm["reversibilite"]  = True

    prm["default"]        = {"tags":"chap","values":["Exercice3_1","py"]} # cas par défaut
    prm["bugs"]           =""" Bugs ou problèmes en cours :
                           aucun apparemment"""
    return prm
 
 
class GuiWin(object):
    """
    Dessin de la fenêtre principale
    """
   
    
    def __init__(self, prm, **kwargs):
        """
        Initial values of the parameters
        """

        self.aide_texte="""
        1 - Sélectionner l'exercice dans la fenêtre principale 
        2 - Si la correction est en FORTRAN, choisir les chemins des dossiers et l'option
        3 - Sélectionner le bouton 'Solution' 
        4 - On peut sélectionner un autre exercice 
        5 - Quitter en appuyant sur 'Sortir' """

        print('#',60*'*')
        print('# %s'%("GUI"))
        print('#',60*'*','\n')
        self.SrcDir         = prm["SrcDir"]             # chemin pour les sources fortran
        self.RunDir         = prm["RunDir"]             # chemin pour lancer la solution en fortran
        self.RunDirPy       = prm["PyRun"]              # Chemin pour sauvegarder les solutions sous python.
     
        self.ShowMenu       = prm["ShowMenu"]           # Utilisation d'un menu
        self.ShowTitre      = prm["ShowTitre"]          # Affichage du titre du livre
        self.ShowTree       = prm["ShowTree"]           # Utilisation d'un tableau
        self.ShowLang       = prm["ShowLang"]           # inutile
        self.ShowFort       = prm["ShowFort"]           # Utilisation du fortran
        self.ShowCheck      = prm["ShowCheck"]          # Utilisation de case  à cocher
        self.ShowText       = prm["ShowText"]           # Position d'un texte (aide) dans la fenêtre principale
        self.ShowDefaultCase = prm["ShowDefaultCase"]    # Affichage du cas par défaut
       
        self.geometry       = prm["geometry"]           # taille de la fenêtre principale
        self.reversibilite  = prm["reversibilite"]
        self.redirection    = prm["redirection"]
        self.original_stdout=sys.stdout
        self.lang           = "Py"                      # "Py" or "F90"
        self.bugs           = prm["bugs"]               # texte pour indiquer les bugs en cours
        #print("parameters : ",prm)

        self.selected_item = False
        self.default       = prm["default"] # cas par défaut
        self.s             = self.default

        if self.s : self.selected_item=True

        self.liste_fortran=liste_exo_fortran(display=False)
        self.liste_python=liste_exo_python(display=False)
        self.liste_chap=liste_chapters()

        self.root = Tk()
        self.root.title("Aero_Suite")
        self.screen_width = self.root.winfo_screenwidth()
        self.screen_height = self.root.winfo_screenheight()
        print("Screen resolution : %i x %i"%(self.screen_width,self.screen_height))
 
 
        if self.ShowMenu: self.def_menu()
        
        #self.root.geometry(self.geometry)   # taille de la fenêtre, optionnel
        self.ftree=Frame(self.root) 
        self.f=Frame(self.root) 

        if self.ShowTree: self.def_tree()
        if self.ShowLang: self.def_Langage()            
        if self.ShowFort: self.def_Fortran()
       
        if self.ShowText: self.def_aide()
        if self.ShowTitre: self.def_Titre()
        if self.ShowDefaultCase : self.def_default()
        if self.ShowCheck: self.def_Check()

        self.stdout2file = IntVar()
        if self.redirection:
            self.stdout2file.set(1)   # ON
        else:
            self.stdout2file.set(0)   # OFF
        Checkbutton(self.f,text="Exo. PYTHON: Sortie standard dans un fichier",variable=self.stdout2file,
            font=("Times","12","bold" ),command=self.change_stdout).grid(row=14,column=0,sticky=W)
        
        self.reverse = IntVar()
        if self.reversibilite:
            self.reverse.set(1)   # ON
        else:
            self.reverse.set(0)   # OFF
        Checkbutton(self.f,text="Reversibilité pour la sortie standard",variable=self.reverse,
            font=("Times","12","bold" ),command=self.change_reverse).grid(row=15,column=0,sticky=W)
   



        # bouton de sortie du code   
        Button(self.f,text="Sortir",bg="tomato",font=("Times","12","bold" ),relief=RAISED,command=self.Sortir, borderwidth=1).grid(row=10,column=0) 
        # bouton pour lancer le calcul
        Button(self.f,text="Calculer Solution",bg="yellow",font=("Times","12","bold" ),relief=SUNKEN,command=self.run).grid(row=11,column=0)


        self.f.pack(side=LEFT)

        # self.fRight=Frame(self.root)

        # # bouton de sortie du code   
        # Button(self.fRight,text="Sortir",bg="tomato",font=("Times","12","bold" ),relief=RAISED,command=self.Sortir, borderwidth=1).grid(row=0,column=4) 
        # # bouton pour lancer le calcul
        # Button(self.fRight,text="Calculer Solution",bg="plum1",font=("Times","12","bold" ),relief=SUNKEN,command=self.run).grid(row=2,column=4)
        # self.fRight.pack(side=BOTTOM)
        self.root.mainloop()

    # ************************************************************
    # FUNCTIONS : 
    # ************************************************************

     
    def def_default(self):
        """
        Affichage du cas par défaut dans la fenêtre principale
        """
        txt_default="Cas par défaut: lang= %s, exo = %s "%(self.s["values"][1],self.s["values"][0])
        Label(self.f, text=txt_default,bg='grey90',font=("Times","11","italic" )).grid(row=13, column=0)

    def fileno(self,file_or_fd):
        """
        Obtenir le numéro du fichier lors de son ouverture
        """
        fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
        if not isinstance(fd, int):
            raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
        return fd

    @contextmanager
    def stdout_redirected(self,to=os.devnull, stdout=None):
        """
        Procédure pour rediriger la sortie standard, l'utilisation du dup2 
        rend la procédure irréversible : on ne peut plus ensuite revenir
        sur la sortie standard habituelle. Il faut relancer le code.
        """
        if stdout is None:
           stdout = sys.stdout

        stdout_fd = self.fileno(stdout)
        # copy stdout_fd before it is overwritten
        #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
        with os.fdopen(os.dup(stdout_fd), 'wb') as copied: 
            stdout.flush()  # flush library buffers that dup2 knows nothing about
            try:
                os.dup2(self.fileno(to), stdout_fd)  # $ exec >&to
            except ValueError:  # filename
                with open(to, 'wb') as to_file:
                    os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
            try:
                yield stdout # allow code to be run with the redirected stdout
            finally:
                # restore stdout to its previous value
                #NOTE: dup2 makes stdout_fd inheritable unconditionally
                stdout.flush()
                os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied
        
    def run(self):
        """
        Opération pour lancer les codes
        """
        sys.stdout = self.original_stdout
        print("redirection : ",self.redirection)
        if self.selected_item:
            print("Résolution d'un exercice")
        else:
            print("aucun élément sélectionné")
            return
#   self.reversibilite=False
# dans un cas, reversibilite, on peut passer modifier la sortie standard plusieurs fois,
#             mais les fichiers sont écrits uniquement après l'opération run suivante 
#             il faut en faire une pour rien

# dans le second cas, les fichiers sont immédiatement remplis, mais on ne peut plus
#                   revenir sur la sortie standard écran. Il faut quitter et revenir.

        if self.s["tags"]!=["chap"]:
            if self.s["values"][1]=="py":
                self.filename=self.dirRunPy.get()+self.s["values"][0]+".out"
                if self.redirection : 
                    self.original_stdout= sys.stdout
                    print("> Nom du fichier : ",self.filename)
                    if self.reversibilite:
                        sys.stdout = open(self.filename, 'w')
                        eval(self.s["values"][0]+'()')
                        # astuce pour remplir le buffer du fichier et le fermer. 
                        # il suffit d'en ouvrir un autre ...
                        sys.stdout = open(str(Path.home())+"/Aeffacer.txt", 'w')
                        print("Fichier à effacer")
                    else:
                        stdout_fd = sys.stdout.fileno()
                        with open(self.filename, 'w') as f, self.stdout_redirected(f):
                            print('redirected to a file')
                            os.write(stdout_fd, b'it is redirected now\n')
                            os.system('echo this is also redirected')
                            try:
                                eval(self.s["values"][0]+'()')
                            except :
                                print('Exception : pas de lancement de ',self.s["values"][0])
                            stdout.flush()
                        sys.stdout = self.original_stdout
                        print("fin de la redirection ?")

                else:
                    sys.stdout = self.original_stdout
                    print("Lancement : ", self.s["values"][0] )
                    eval(self.s["values"][0]+'()')
                 
            
            elif self.s["values"][1]=="f90":
                if platform.system()=="win32":
                    self.root.option_add("*Dialog.msg.width",35)
                    txt="""
                    Le lancement 
                    automatique 
                    du code FORTRAN
                    n'est pas implémenté 
                    pour WINDOWS
                    """
                    showwarning("Problème",txt)
                else:
                    print("code en fortran");
                    for k,liste in enumerate(liste_exo_fortran(display=False)):
                        if self.s["values"][0] == liste['exo']:
                            print("valeur trouvée = ", k)
                            print("index de la liste envoyé :", [k])
                            run_fortran(self.opt_fortran.get(),[k],SrcDir=self.dirSrc.get(),RunDir=self.dirRun.get())
                    
            #exit()
        # il semblerait qui n'arrive pas ici
         
        # if self.redirection:
        sys.stdout = self.original_stdout
        #     print("end run",sys.stdout)
        print("Fin du calcul : nouveau calcul ?")    
    
    def Show_Solution(self,Lang):
        """
        Afficher la solution dans une fenêtre
        pour le fortran uniquement
        """
        def Quit():
            fen.destroy()

        if self.redirection:
            # astuce pour remplir le buffer du fichier et le fermer. 
            # il suffit d'en ouvrir un autre ...
            sys.stdout = open(str(Path.home())+"/Aeffacer.txt", 'w')
            print("Fichier à effacer")
                

        fen = Tk()
        if Lang=='Py':
            self.filename=filedialog.askopenfilename(initialdir=self.dirRunPy.get(),
                filetypes=[("output files","*.out"),("data files","*.dat"),('text files', '.txt'), ('all files', '.*')])
        else:
            self.filename=filedialog.askopenfilename(initialdir=self.dirRun.get(),
                filetypes=[("output files","*.out"),("data files","*.dat"),('text files', '.txt'), ('all files', '.*')])
        fen.geometry("600x500+300+300")
        fen.grid_rowconfigure(0,weight=1)
        fen.grid_columnconfigure(0,weight=1)
        fen.title("Fichier")

        menubar=Menu(fen)
        filemenu=Menu(menubar,tearoff=0)
        filemenu.add_command(label="Sortir",command=Quit)
        menubar.add_cascade(label="Fichier",menu=filemenu)
        fen.config(menu=menubar)        
                    
        hscrollbar=Scrollbar(fen,orient=HORIZONTAL)
        hscrollbar.grid(row=1,column=0,sticky=E+W)

        vscrollbar=Scrollbar(fen)
        vscrollbar.grid(row=0,column=1,sticky=N+S)

         
        text=Text(fen,wrap=NONE,bg='light yellow',bd=1,xscrollcommand=hscrollbar.set,yscrollcommand=vscrollbar.set)
        text.grid(row=0,column=0,sticky=E+W+S+N)

        text.insert("1.0","Choisir le fichier dans le menu")
        hscrollbar.config(command=text.xview)
        vscrollbar.config(command=text.yview)

        if self.filename:
            fen.title(self.filename)
            #original=text.get("0.0","end-1c")
            with open(self.filename) as filetoread:
                txtfilecontent=filetoread.read()
                filetoread.close()
            text.delete("1.0", END) # Erase the previous Text widget content
            text.insert("1.0", txtfilecontent)
       
        fen.mainloop()

    def change_reverse(self):
        if self.reverse.get()==1:
            self.reversibilite=True
        else:
            self.reversibilite=False

    def change_stdout(self):
        if self.stdout2file.get()==1:
            self.redirection=True
            print("redirection des sorties")
        else:
            self.redirection=False
            sys.stdout = self.original_stdout
            print("sortie standard")  
        print("fin change, redirection :",self.redirection)          

   
    def def_Fortran(self):
        """
        Bloc pour lancer des exercices écrits en fortran
        """
        #Options pour le fortran  
        options_fortran=["compilation","run","compil+run"]
        Label(self.f, text="Option du Fortran : ",font=("Times","11","bold" )).grid(row=0,column=0,sticky=W)
        self.opt_fortran=StringVar()
        self.opt_fortran.set("run")
        OptionMenu(self.f, self.opt_fortran,*options_fortran).grid(row=1,column=0,sticky=W)

        # chemins par défaut
        # pour les sources Fortran et le dossier de l'éxécutable
        self.dirSrc=StringVar();  self.dirSrc.set(self.SrcDir)
        self.dirRun=StringVar();  self.dirRun.set(self.RunDir)
        self.dirRunPy=StringVar();self.dirRunPy.set(self.RunDirPy)
        print("dirRun",self.dirRun.get())
        # Boutons pour sortie, lancer le code, choisir les dossiers Fortran

        # entrée des dossiers FORTRAN
        Button(self.f,textvar=self.dirSrc,bg="plum1", command=self.browse_source).grid(row=3, column=0,sticky=W)
        Label(self.f, text="Fortran, dossier SOURCE : ",bg='white',font=("Times","11","bold" )).grid(row=2, column=0,sticky=W)
        Button(self.f,textvar=self.dirRun,bg="plum1", command=self.browse_run).grid(row=5, column=0,sticky=W)
        Label(self.f, text="Fortran, dossier RUN : ",bg='white', font=("Times","11","bold" )).grid(row=4, column=0,sticky=W)
        Button(self.f,textvar=self.dirRunPy,bg="skyblue", command=self.browse_runPy).grid(row=7, column=0,sticky=W)
        Label(self.f, text="Python, dossier RUN : ",bg='white', font=("Times","11","bold" )).grid(row=6, column=0,sticky=W)


    def def_tree(self):

        """
        tableau des exercices, 
        """
        self.tree = Treeview(self.root,height=17,selectmode="browse")

        style = ttk.Style(self.root)
        style.configure('Treeview', rowheight=30, background="grey90", foreground="black", fieldbackground="grey99")

        vsb = ttk.Scrollbar(self.tree,orient="vertical",command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        vsb.place(x=745, y=20, height=500)

        #tree.configure(yscrollcommand=vsb.set)
        self.tree["columns"]=("one","two","three")
        self.tree.column("#0", width=100) 
        self.tree.column("one", width=110 )
        self.tree.column("two", width=50)
        self.tree.column("three", width=500)
        self.tree.heading("#0",text="Livre")
        self.tree.heading("one", text="Numéro")
        self.tree.heading("two", text="Lang")
        self.tree.heading("three", text="Commentaires")
        k=0
        ch_list=[]
        name_list=[]
        for i in range(15):
            ch_list.append("ch. " + str(1+i))
            name_list.append("dir" + str(1+i))
        
          
        for txt,name,chap in zip(ch_list,name_list,self.liste_chap):
            self.tree.insert("", k, name, text=txt,tags="chap",values=("","",chap))
            #tree.insert(name, k, text="Exercice",values=("3A"," 3B"))
            #tree.insert(name, k, text="Problème",values=("3A"," 3B"))
            #for j in range(3):
            #     tree.insert(name, k, text="Exo",values=("12-"+str(j)," 3B"))
            k=k+1
            
        for liste in self.liste_fortran:
            k=liste["ch"]-1
            self.tree.insert(name_list[k],k,text="Exo",values=(liste['exo'],"f90",liste['com']),tags="exof90")
        for liste in self.liste_python:
            k=liste["ch"]-1
            self.tree.insert(name_list[k],k,text="Exo",values=(liste['exo'],"py",liste['com']),tags="exopy")
        self.tree.tag_configure('exof90', background='plum1')
        self.tree.tag_configure('exopy', background='LightSkyBlue1')
        #tree.grid(row=10, column=0, sticky='nsew')
        self.tree.bind('<<TreeviewSelect>>', self.on_selected)
        
        self.tree.pack(side=LEFT)


    def run_default_case(self):
        """
        pour lancer la cas par défaut
        """
        print("run the default case") 
        self.s=self.default
        print("s = ",self.s)
        self.selected_item=True



    def def_menu(self):
        """
        DEFINITION DE MENUS DE LA FENETRE
        """
        menubar = Menu(self.root,bg='grey50')

        menu1 = Menu(menubar, tearoff=0, background='#525049', foreground='white',
               activebackground='white', activeforeground='black')
        menu1.add_command(label="cas par défaut", font=("Ubuntu", 10),command=self.run_default_case)
        menu1.add_command(label="A faire & bugs", font=("Ubuntu", 10),command=self.Show_Bugs)
        menu1.add_separator()
        menu1.add_command(label="Quitter",font=("Ubuntu", 10), command=self.root.quit)
        menubar.add_cascade(label="Général",font=("Ubuntu", 10), menu=menu1)

        menu2 = Menu(menubar, tearoff=0,background='#525049', foreground='white',
               activebackground='white', activeforeground='black')
        menu2.add_command(label="Affichage",font=("Ubuntu", 10), command=self.alert)
        menu2.add_command(label="Fortran", font=("Ubuntu", 10),command=self.alert)
        menu2.add_command(label="Python", font=("Ubuntu", 10),command=self.alert)
        menubar.add_cascade(label="Options", font=("Ubuntu", 10),menu=menu2)

        menu3 = Menu(menubar, tearoff=0,background='#525049', foreground='white',
               activebackground='white', activeforeground='black')
        menu3.add_command(label="Calculer la solution",font=("Ubuntu", 10), command=self.run)
        menu3.add_command(label="Montrer Solution Fortran", font=("Ubuntu", 10),command=lambda:self.Show_Solution("F90"))
        menu3.add_command(label="Montrer Solution Python",font=("Ubuntu", 10), command=lambda: self.Show_Solution("Py"))
        menu3.add_command(label="Effacer Solution Python", font=("Ubuntu", 10),command=self.alert)
        menubar.add_cascade(label="Solution", font=("Ubuntu", 10),menu=menu3)


        menu4 = Menu(menubar, tearoff=0,background='#525049', foreground='white',
               activebackground='white', activeforeground='black')
        menu4.add_command(label="A propos", font=("Ubuntu", 10),command=self.apropos)
        menu4.add_command(label="Aide", font=("Ubuntu", 10),command=self.aide)
        menubar.add_cascade(label="Aide",font=("Ubuntu", 10), menu=menu4)

        self.root.config(menu=menubar)



    def Show_Bugs(self):
        """
        Informations sur les bugs, 
        ou les choses à faire.
        """
        if len(self.bugs)>=0:
            self.display_info()
        else: 
            showinfo("Information bugs / à faire ", "rien")

    def display_info(self):
        """
        Message d'information dans une fenêtre
        """
        def Quit():
            f_info.destroy()
        
        """
        Menu Informations bugs
        """
        f_info =  Tk()
        f_info.geometry('640x300+40+40')   # taille de la fenêtre, optionnel
        F=Frame(f_info)
        f_info.title("Information 'à faire' ou 'bugs'")
        T=Text(f_info,height=6, width=100,font=("Times","14","italic" ))
        T.pack()
        T.insert(END,self.bugs)
        Button(f_info,text="Sortir",bg="tomato",font=("Times","12","bold" ),command=Quit).pack(side="bottom")
        F.pack()
        f_info.mainloop()


    def def_Titre(self):
        """
        Affichage du titre du livre
        """
        Titre=Text(self.root,bg='bisque',height=6, width=80)
        Titre.pack()
        txt="""Suite Aérodynamique : 
        solution des exercices et problèmes du livre
        'Exercices et Problèmes d'Aérodynamique',
        auteur: Christophe Airiau, André Giovannini, Pierre Brancher
        Université Paul Sabatier, Toulouse III
        Institut de Mécanique des Fluides de Toulouse.
        """
        Titre.insert(END,txt)



    def def_aide(self):
        """affichage d'un message """
        #Label(f,text='je place le texte où je veux').place(x=-10,y=50)
        #Label(f,text='je place le texte où je veux').grid(row=4,column=0)
        T=Text(self.root,height=6, width=140)
        T.pack()
        T.insert(END,self.aide_texte)

    def def_Check(self):
        """
        boite check d'option : affichage de liste
        """
        self.f.pack()
        liste_exos_python = IntVar()
        liste_exos_python.set(1)   # ON
        Checkbutton(self.f,text="Liste exos Python",variable=liste_exos_python,command=self.liste_python).grid(row=14,column=0)
        liste_exos_fortran = IntVar()
        liste_exos_fortran.set(0)   # ON
        Checkbutton(self.f,text="Liste exos Fortran",variable=liste_exos_fortran,command=self.liste_fortran).grid(row=15,column=0)


    def def_Langage(self):
        """
        Pour montrer un exemple de case à cocher
        mais pas utiliser en l'état
        """        
 
        choix=["Python","Fortran"]
        var=StringVar()
        var.set("Python")
        Label(self.f, text="Choix du langage : ",font=("Times","11","bold" )).grid(row=17,column=0)
        OptionMenu(self.f,var,*choix).grid(row=7,column=1)


    def aide(self):
        """
        Affiche de l'aide au programme
        """
        def Quit():
            f_aide.destroy()
        
        """
        Menu Aide
        """
        f_aide = Tk()
        f_aide.geometry('700x200+40+40')   # taille de la fenêtre, optionnel
        F=Frame(f_aide)
        f_aide.title("Aide")
        T=Text(f_aide,height=6, width=140,font=("Times","12","italic" ))
        T.pack()
        T.insert(END,self.aide_texte)
        Button(f_aide,text="Sortir",bg="tomato",font=("Times","12","bold" ),command=Quit).pack(side="bottom")
        F.pack()
        f_aide.mainloop()

    def apropos(self):
        """ A propos box """
        txt=""" Suite Aérodynamique : 
        Solutions trouvées dans le  livre 
        'Exercices d'Aérodynamique Fondamentale',
        C. Airiau,
        A. Giovannini,
        P. Brancher
        Codes écrit par C. Airiau
        """
        self.root.option_add("*Dialog.msg.width",35)
        showinfo("A propos", txt)
            

 

    def on_selected(self,x):
        """
        details of the selected item
        """
        self.selected_item=True
        self.s=self.tree.item(self.tree.focus())
        # print("selection : ",s)
         
        # if s["tags"]!=["chap"] :
        #     print("codage  : ",s["tags"])
        #     print("commentaire : ",s["values"][2])
        #     if s["values"][1]=="f90":
        #         for k,liste in enumerate(liste_exo_fortran(display=False)):
        #             if s["values"][0] == liste['exo']:
        #                 print("valeur trouvée = ", k)

    def alert(self):
        """ Alert box """
        showinfo("alerte", "A implémenter")

    def browse_source(self):
        self.dirSrc.set(filedialog.askdirectory())
    
    def browse_run(self):
        self.dirRun.set(filedialog.askdirectory())

    def browse_runPy(self):
        self.dirRunPy.set(filedialog.askdirectory())
     
    def Sortir(self):
        print("Quitter")
        try:
            os.remove(str(Path.home())+"/Aeffacer.txt")
        except:
            pass
        finally:
            exit() 

    def liste_python():
        """
        utiliser en cas de liste à cocher, inutile pour le moment
        """
        print("afficher la liste des programmes en python")

    def liste_fortran():
        """
        utiliser en cas de liste à cocher, inutile pour le moment
        """
        print("afficher la liste des programmes en fortran")
