#!/bin/py
# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2020/04/10*

*licence : AGPL-3.0*

Module to get the solutions of the exercises using a tkinter GUI 
  ..
"""
from contextlib import contextmanager
import os
from sys import stdout
import platform
# import tkinter
from tkinter import Tk, Frame, Checkbutton, IntVar, StringVar, Button, PhotoImage, Canvas
from tkinter import RAISED, SUNKEN, LEFT, Label, Menu, Scrollbar, HORIZONTAL
from tkinter import E, N, S, W, Text, NONE, END, OptionMenu, TOP
from tkinter.messagebox import showwarning, showinfo
from tkinter import ttk
from tkinter import filedialog


from default_path import PyRun, FortranRun, FortranPath
from Book.Fortran.fortran import exercise_list_fortran, run_EC,  run_fortran, sys, Path
import Tools.misc

from CompressibleFlow.fonctions import *
from Book.Chapter2 import *
from Book.Chapter3 import *
from Book.Chapter4 import *
from Book.Chapter5 import *
from Book.Chapter6 import *
from Book.Chapter7 import *
from Book.Chapter7c import *
from Book.Chapter8 import *
from Book.Chapter9 import *
from Book.Chapter10 import *
from Book.Chapter11 import *
from Book.Chapter12 import *
from Book.Chapter12b import *
from Book.Chapter13 import *
from Book.Chapter14 import *
from Book.Chapter15_Laminaire import Exercice15_1, Exercice15_2, Exercice15_4, Exercice15_6
from Book.Chapter15_Laminaire import Exercice15_7, Exercice15_8
from Book.Chapter15_Transition import Exercice15_9
from Book.Chapter15_Turbulent import Exercice15_10
from Book.Chapter15_Tu import Exercice15_11, Exercice15_12, Exercice15_13
from Book.Chapter15_Tu import Exercice15_14, Exercice15_15
from Book.Corrections import liste_exo_python, liste_chapters
from Book.misc import *
from conical_shock.conical_shock import *
from Book.panneaux_intro import *

def set_params_gui():
    """
    Set default parameters
        ..
    """
    prm = {}

    prm["SrcDir"] = FortranPath             # path to the Fortran source file
    prm["RunDir"] = FortranRun              # Path to run  the solution with Fortran
    prm["PyRun"] = PyRun                    # Path to save the solution with Python in an output file
    prm["ShowMenu"] = True                  # to use a menu
    prm["ShowTree"] = True                  # to use a table
    prm["ShowLang"] = False                 # unused ?
    prm["ShowFort"] = True                  # to use Fortran
    prm["ShowCheck"] = False                # to use check square
    prm["ShowText"] = False                 #
    prm["ShowTitre"] = False                # to diplay the book title
    prm["ShowDefaultCase"] = True           # location of the help text in the main window
    prm["geometry"] = '1000x800+100+50'     # size of the main window
    prm["redirection"] = False
    prm["reversibility"] = True
    prm["language"] = "en"                  # "en" or "fr" ("fr" not updated)
    prm["default"] = {"tags": "chap", "values": ["Exercice3_1", "py"]}  # default case
    prm["bugs"] = """ Bugs or current problem  : 
        Depending of how you use the GUI (spyder, into a terminal, pycharm, etc ...) 
        the behaviour of the exit or quit menu works well or not,
        or it is possible to run several or only one exercise solution.
        
        I have some understanding problem with TkInter on the window management.
        
        The solutions have been validated each time it was possible.
        But unfortunately I can not warranty that no bug remains ! 
    """

    print("Parameters : default fortran run dir : ", prm["RunDir"])
    return prm


class GuiWin(object):
    """
    class to generate the main window with the menu
    """

    def __init__(self, prm, **kwargs):
        """
        Initial values of the parameters
        """
        if prm["language"] == "fr":
            self.aide_texte = """
            1 - Sélectionner l'exercice dans la fenêtre principale 
            2 - Si la correction est en FORTRAN, choisir les chemins des dossiers et l'option
            3 - Sélectionner le bouton 'Solve' 
            4 - On peut sélectionner un autre exercice 
            5 - Quitter en appuyant sur 'Exit' """
        else:
            self.aide_texte = """
            1 - Select the exercise in the main window 
            2 - If the solution is with Fortran please enter the Fortran directory paths and the option
            3 - Select the 'Solve' button
            4 - Then one can select another exercise 
            5 - Quit with the 'Exit' button """

        print('#', 60 * '*')
        print('# %s' % ("GUI"))
        print('#', 60 * '*', '\n')
        self.SrcDir = prm["SrcDir"]                 # Path for Fortran sources
        self.RunDir = prm["RunDir"]                 # Path to run  the solution with Fortran
        self.RunDirPy = prm["PyRun"]                # Path to save the solution with Python in an output file
        self.ShowMenu = prm["ShowMenu"]             # to use a menu
        self.ShowTitre = prm["ShowTitre"]           # to diplay the book title
        self.ShowTree = prm["ShowTree"]             # to use a table
        self.ShowLang = prm["ShowLang"]             # unused ?
        self.ShowFort = prm["ShowFort"]             # to use Fortran
        self.ShowCheck = prm["ShowCheck"]           # to use check square
        self.ShowText = prm["ShowText"]             # location of the help text in the main window
        self.ShowDefaultCase = prm["ShowDefaultCase"]  # to display the default case
        self.geometry = prm["geometry"]             # size of the main window
        self.reversibility = prm["reversibility"]   # reversibility in the ouput redirection
        self.redirection = prm["redirection"]       # for redirection of the output from the terminal to a file
        self.original_stdout = sys.stdout           # standard output
        self.lang = "Py"                            # "Py" or "F90"
        self.bugs = prm["bugs"]                     # text for the current unsolved bugs
        self.language = prm["language"]             # "en" or "fr" (just for the title of exercises)
        self.show_logo = True

        self.selected_item = False
        self.default = prm["default"]  # default case
        self.s = self.default
        print("help: ", self.aide_texte)

        if self.s:
            self.selected_item = True
            print("defaut exercise selected : ", self.s["values"])

        self.liste_fortran = exercise_list_fortran(display=False)
        self.liste_python = liste_exo_python(display=False, lang=self.language)
        self.liste_chap = liste_chapters(lang=self.language)

        self.root = Tk()
        self.root.title("FundAeroSuite")
        self.screen_width = self.root.winfo_screenwidth()
        self.screen_height = self.root.winfo_screenheight()
        print("Screen resolution : %i x %i" % (self.screen_width, self.screen_height))

        print("bugs : ", self.bugs)
        if self.ShowMenu:
            self.def_menu()

        # self.root.geometry(self.geometry)             # size of the main windows, optional
        self.ftree = Frame(self.root)
        self.f = Frame(self.root)

        if self.ShowTree:
            self.def_tree()
        if self.ShowLang:
            self.def_Langage()
        if self.ShowFort:
            self.def_Fortran()
        if self.ShowText:
            self.def_aide()
        if self.ShowTitre:
            self.def_Titre()
        if self.ShowDefaultCase:
            self.def_default()
        if self.ShowCheck:
            self.def_Check()

        self.stdout2file = IntVar()
        if self.redirection:
            self.stdout2file.set(1)  # ON
        else:
            self.stdout2file.set(0)  # OFF
        Checkbutton(self.f, text="Exo. PYTHON: standard output in a file", variable=self.stdout2file,
                    font=("Times", "12", "bold"), command=self.change_stdout).grid(row=14, column=0, sticky=W)

        self.reverse = IntVar()
        if self.reversibility:
            self.reverse.set(1)  # ON
        else:
            self.reverse.set(0)  # OFF
        Checkbutton(self.f, text="Reversibility of the standard output", variable=self.reverse,
                    font=("Times", "12", "bold"), command=self.change_reverse).grid(row=15, column=0, sticky=W)

        # Exit button
        Button(self.f, text="Exit <Ctrl-q>", bg="tomato", font=("Times", "12", "bold"), relief=RAISED, command=self.Sortir,
               borderwidth=1).grid(row=10, column=0)
        # button to get solution of an exercise
        Button(self.f, text="Solve <Ctrl-s>", bg="yellow", font=("Times", "12", "bold"), relief=SUNKEN, command=self.run).grid(
            row=11, column=0)

        if self.show_logo :
            photo = PhotoImage(file="image/logo_doc.png")
            label = Label(self.f, image=photo, width=320, height=123, relief=SUNKEN).grid(row=1, column=1)
        self.f.pack(side=LEFT)
            # canvas = Canvas(self.f, width=128, height=128, bg="tomato").grid(row=10, column=1)

        self.root.bind('<Control-q>',self.Sortir)
        self.root.bind('<Control-s>', self.run)
        self.root.mainloop()
        print("end of GUI execution")

    # ************************************************************
    # FUNCTIONS : 
    # ************************************************************

    def def_default(self):
        """
        Display the default case in the main window
        """
        txt_default = "default case: lang= %s, exo = %s " % (self.s["values"][1], self.s["values"][0])
        Label(self.f, text=txt_default, bg='grey90', font=("Times", "11", "italic")).grid(row=13, column=0)

    def fileno(self, file_or_fd):
        """
        get the file number when opening
        """
        fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
        if not isinstance(fd, int):
            raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
        return fd

    @contextmanager
    def stdout_redirected(self, to=os.devnull, stdout=None):
        """
        procedure to redirect the standard output
        The use of the dup2 render the procedure irreversible. We can't no longer come
        back to the usual standard output (terminal).
        We have to open the main script again 
        """
        if stdout is None:
            stdout = sys.stdout

        stdout_fd = self.fileno(stdout)
        # copy stdout_fd before it is overwritten
        # NOTE: `copied` is inheritable on Windows when duplicating a standard stream
        with os.fdopen(os.dup(stdout_fd), 'wb') as copied:
            stdout.flush()  # flush library buffers that dup2 knows nothing about
            try:
                os.dup2(self.fileno(to), stdout_fd)  # $ exec >&to
            except ValueError:  # filename
                with open(to, 'wb') as to_file:
                    os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
            try:
                yield stdout  # allow code to be run with the redirected stdout
            finally:
                # restore stdout to its previous value
                # NOTE: dup2 makes stdout_fd inheritable unconditionally
                stdout.flush()
                os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

    def run(self, event=""):
        """
        function to run the main code
        """
        sys.stdout = self.original_stdout
        print("redirection : ", self.redirection)
        if self.selected_item:
            print("Solve an exercise")
        else:
            print("no element is selected")
            return
        #   self.reversibility=False
        #
        # dans un cas, reversibility, on peut passer modifier la sortie standard plusieurs fois,
        #             mais les fichiers sont écrits uniquement après l'opération run suivante
        #             il faut en faire une pour rien
        # dans le second cas, les fichiers sont immédiatement remplis, mais on ne peut plus
        #                   revenir sur la sortie standard écran. Il faut quitter et revenir.

        print('selection : ', self.selected_item)
        if self.s["tags"] != ["chap"]:
            if self.s["values"][1] == "py":
                self.filename = self.dirRunPy.get() + "/" + self.s["values"][0] + ".out"
                if self.redirection:
                    self.original_stdout = sys.stdout
                    print("> File name : ", self.filename)
                    if self.reversibility:
                        sys.stdout = open(self.filename, 'w')
                        eval(self.s["values"][0] + '()')
                        #  a way to fill the file buffer and to close it. 
                        # we have to open a new one ...
                        sys.stdout = open(str(Path.home()) + "/ToDelete.txt", 'w')
                        print("File to delete")
                    else:
                        stdout_fd = sys.stdout.fileno()
                        with open(self.filename, 'w') as f, self.stdout_redirected(f):
                            print('redirected to a file')
                            os.write(stdout_fd, b'it is redirected now\n')
                            os.system('echo this is also redirected')
                            try:
                                eval(self.s["values"][0] + '()')
                            except:
                                print('Exception : no run of  ', self.s["values"][0])
                            stdout.flush()
                        sys.stdout = self.original_stdout
                        print("redirection end ?")

                else:
                    sys.stdout = self.original_stdout
                    print("run : ", self.s["values"][0])
                    eval(self.s["values"][0] + '()')

            elif self.s["values"][1] == "f90":
                if platform.system() == "win32":
                    self.root.option_add("*Dialog.msg.width", 35)
                    txt = """
                    to run automatically the Fortran sources is not implemented in WINDOWS
                    """
                    showwarning("Problem : ", txt)
                else:
                    print("*** Fortran code ***")
                    print('selected exercise :', self.s['values'][0])
                    for k, liste in enumerate(exercise_list_fortran(display=False)):
                        if self.s["values"][0] == liste['exo']:
                            print("found value  = ", k," for chap = ", liste["ch"], " n° ", liste["num"] )
                            exo = [[liste["ch"], [liste["num"]]]]
                            run_fortran(self.opt_fortran.get(), exo, SrcDir=self.dirSrc.get(), RunDir=self.dirRun.get())

            # exit()
        # il semblerait qu'on n'arrive jamais ici: est-ce un problème ?

        # if self.redirection:
        sys.stdout = self.original_stdout
        #     print("end run",sys.stdout)
        print("End of the calculus : to try another one use the selector")

    def Show_Solution(self, Lang):
        """
        To display the solution into a new window,
            * mainly for  Fortran use
            * possibly for Python solution if the outputs have been redirected into files
        """

        def Quit():
            fen.destroy()

        if self.redirection:
            #  a way to fill the file buffer and to close it. 
            # we have to open a new one ...
            sys.stdout = open(str(Path.home()) + "/ToDelete.txt", 'w')
            print("File to delete")

        fen = Tk()
        if Lang == 'Py':
            self.filename = filedialog.askopenfilename(initialdir=self.dirRunPy.get(),
                                                       filetypes=[("output files", "*.out"), ("data files", "*.dat"),
                                                                  ('text files', '.txt'), ('all files', '.*')])
        else:
            self.filename = filedialog.askopenfilename(initialdir=self.dirRun.get(),
                                                       filetypes=[("output files", "*.out"), ("data files", "*.dat"),
                                                                  ('text files', '.txt'), ('all files', '.*')])
        fen.geometry("600x500+300+300")
        fen.grid_rowconfigure(0, weight=1)
        fen.grid_columnconfigure(0, weight=1)
        fen.title("File")

        menubar = Menu(fen)
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="Exit", command=Quit)
        menubar.add_cascade(label="File", menu=filemenu)
        fen.config(menu=menubar)

        hscrollbar = Scrollbar(fen, orient=HORIZONTAL)
        hscrollbar.grid(row=1, column=0, sticky=E + W)
        vscrollbar = Scrollbar(fen)
        vscrollbar.grid(row=0, column=1, sticky=N + S)

        text = Text(fen, wrap=NONE, bg='light yellow', bd=1, xscrollcommand=hscrollbar.set,
                    yscrollcommand=vscrollbar.set)
        text.grid(row=0, column=0, sticky=E + W + S + N)

        text.insert("2.0", "Choose a file with the menu")
        hscrollbar.config(command=text.xview)
        vscrollbar.config(command=text.yview)

        if self.filename:
            fen.title(self.filename)
            # original=text.get("0.0","end-1c")
            with open(self.filename) as filetoread:
                txtfilecontent = filetoread.read()
                filetoread.close()
            text.delete("1.0", END)  # Erase the previous Text widget content
            text.insert("1.0", txtfilecontent)

        fen.mainloop()

    def change_reverse(self):
        if self.reverse.get() == 1:
            self.reversibility = True
        else:
            self.reversibility = False

    def change_stdout(self):
        if self.stdout2file.get() == 1:
            self.redirection = True
            print("Outputs redirection")
        else:
            self.redirection = False
            sys.stdout = self.original_stdout
            print("Standard outut")
        print("end of change, redirection :", self.redirection)

    def def_Fortran(self):
        """
        function to run the exercise with Fortran
        """
        # Fortran options
        options_fortran = ["compilation", "run", "compil+run"]
        Label(self.f, text="Fortran option : ", font=("Times", "11", "bold")).grid(row=0, column=0, sticky=W)
        self.opt_fortran = StringVar()
        self.opt_fortran.set("compil+run")
        OptionMenu(self.f, self.opt_fortran, *options_fortran).grid(row=1, column=0, sticky=W)

        # default paths
        # for Fortran source  and the directory of the binary code
        self.dirSrc = StringVar()
        self.dirSrc.set(self.SrcDir)
        self.dirRun = StringVar()
        self.dirRun.set(self.RunDir)
        self.dirRunPy = StringVar()
        self.dirRunPy.set(self.RunDirPy)
        print("dirRun", self.dirRun.get())
        # Output button, run the code, choose the Fortran directory

        #  FORTRAN directories : 
        Button(self.f, textvar=self.dirSrc, bg="plum1", command=self.browse_source).grid(row=3, column=0, sticky=W)
        Label(self.f, text="Fortran, SOURCE directory: ", bg='white', font=("Times", "11", "bold")).grid(row=2,
                                                                                                         column=0,
                                                                                                         sticky=W)
        Button(self.f, textvar=self.dirRun, bg="plum1", command=self.browse_run).grid(row=5, column=0, sticky=W)
        Label(self.f, text="Fortran, RUN directory : ", bg='white', font=("Times", "11", "bold")).grid(row=4, column=0,
                                                                                                       sticky=W)
        Button(self.f, textvar=self.dirRunPy, bg="skyblue", command=self.browse_runPy).grid(row=7, column=0, sticky=W)
        Label(self.f, text="Python, RUN directory : ", bg='white', font=("Times", "11", "bold")).grid(row=6, column=0,
                                                                                                      sticky=W)

    def def_tree(self):
        """
        exercise table 
        """
        self.tree = ttk.Treeview(self.root, height=17, selectmode="browse")
        style = ttk.Style(self.root)
        style.configure('Treeview', rowheight=30, background="grey90", foreground="black", fieldbackground="grey99")
        vsb = ttk.Scrollbar(self.tree, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        vsb.place(x=745, y=20, height=500)

        # tree.configure(yscrollcommand=vsb.set)
        self.tree["columns"] = ("one", "two", "three")
        self.tree.column("#0", width=100)
        self.tree.column("one", width=110)
        self.tree.column("two", width=50)
        self.tree.column("three", width=500)
        self.tree.heading("#0", text="Book")
        self.tree.heading("one", text="Number")
        self.tree.heading("two", text="Lang")
        self.tree.heading("three", text="Comment")
        k = 0
        ch_list = []
        name_list = []
        for i in range(15):
            ch_list.append("ch. " + str(1 + i))
            name_list.append("dir" + str(1 + i))

        for txt, name, chap in zip(ch_list, name_list, self.liste_chap):
            self.tree.insert("", k, name, text=txt, tags="chap", values=("", "", chap))
            k += 1

        for liste in self.liste_fortran:
            k = liste["ch"] - 1
            self.tree.insert(name_list[k], k, text="Exo", values=(liste['exo'], "f90", liste['com']), tags="exof90")
        for liste in self.liste_python:
            k = liste["ch"] - 1
            self.tree.insert(name_list[k], k, text="Exo", values=(liste['exo'], "py", liste['com']), tags="exopy")
        self.tree.tag_configure('exof90', background='plum1')
        self.tree.tag_configure('exopy', background='LightSkyBlue1')
        # tree.grid(row=10, column=0, sticky='nsew')
        self.tree.bind('<<TreeviewSelect>>', self.on_selected)

        self.tree.pack(side=LEFT)

    def run_default_case(self):
        """
        to run the default case
        """
        print("run the default case")
        self.s = self.default
        print("s = ", self.s)
        self.selected_item = True

    def def_menu(self):
        """
        Definition of the menus in the window
        """
        menubar = Menu(self.root, bg='grey50')

        menu1 = Menu(menubar, tearoff=0, background='#525049', foreground='white',
                     activebackground='white', activeforeground='black')
        menu1.add_command(label="default case", font=("Ubuntu", 10), command=self.run_default_case)
        menu1.add_command(label="To do  & bugs", font=("Ubuntu", 10), command=self.Show_Bugs)
        menu1.add_separator()
        menu1.add_command(label="Quit", font=("Ubuntu", 10), command=self.root.quit, accelerator="Ctrl+Q")
        menubar.add_cascade(label="General", font=("Ubuntu", 10), menu=menu1)

        # menu2 = Menu(menubar, tearoff=0, background='#525049', foreground='white',
        #              activebackground='white', activeforeground='black')
        # menu2.add_command(label="Display", font=("Ubuntu", 10), command=self.alert)
        # menu2.add_command(label="Fortran", font=("Ubuntu", 10), command=self.alert)
        # menu2.add_command(label="Python", font=("Ubuntu", 10), command=self.alert)
        # menubar.add_cascade(label="Options", font=("Ubuntu", 10), menu=menu2)

        menu3 = Menu(menubar, tearoff=0, background='#525049', foreground='white',
                     activebackground='white', activeforeground='black')
        menu3.add_command(label="Solve", font=("Ubuntu", 10), command=self.run)
        menu3.add_command(label="Display the Fortran Solution", font=("Ubuntu", 10),
                          command=lambda: self.Show_Solution("F90"))
        menu3.add_command(label="Display the Python Solution", font=("Ubuntu", 10),
                          command=lambda: self.Show_Solution("Py"))
        menu3.add_command(label="Delete the Python Solution", font=("Ubuntu", 10), command=self.alert)
        menubar.add_cascade(label="Solution", font=("Ubuntu", 10), menu=menu3)

        menu4 = Menu(menubar, tearoff=0, background='#525049', foreground='white',
                     activebackground='white', activeforeground='black')
        menu4.add_command(label="About", font=("Ubuntu", 10), command=self.apropos)
        menu4.add_command(label="Help", font=("Ubuntu", 10), command=self.aide)
        menubar.add_cascade(label="Help", font=("Ubuntu", 10), menu=menu4)

        self.root.config(menu=menubar)

    def Show_Bugs(self):
        """
        Informations about the bugs, 
        or the things to do
        """
        if not self.bugs:
            self.display_info()
        else:
            showinfo("Information bugs / to do ", self.bugs)

    def display_info(self):
        """
        Information message in a window
        """

        def Quit():
            f_info.destroy()

        """
        Menu bug Information
        """
        f_info = Tk()
        f_info.geometry('640x300+40+40')  # window size, optional
        F = Frame(f_info)
        f_info.title("'to do' or 'bugs' information")
        T = Text(f_info, height=6, width=100, font=("Times", "14", "italic"))
        T.pack()
        T.insert(END, self.bugs)
        Button(f_info, text="Exit", bg="tomato", font=("Times", "12", "bold"), command=Quit()).pack(side="bottom")
        F.pack()
        f_info.mainloop()

    def def_Titre(self):
        """
        to diplay the book title
        """
        Titre = Text(self.root, bg='bisque', height=6, width=80)
        Titre.pack()
        txt = """Aerodynamic suite : 
        solution of the exercises and problem of the book 
        'Exercices et Problèmes d'Aérodynamique',
        authors: Christophe Airiau, André Giovannini, Pierre Brancher
        Université Paul Sabatier, Toulouse III
        Institut de Mécanique des Fluides de Toulouse.
        """
        Titre.insert(END, txt)

    def def_aide(self):
        """to display messages """
        # Label(f,text='je place le texte où je veux').place(x=-10,y=50)
        # Label(f,text='je place le texte où je veux').grid(row=4,column=0)
        T = Text(self.root, height=6, width=140)
        T.pack()
        T.insert(END, self.aide_texte)

    def def_Check(self):
        """
        check boxes : to display list  
        """
        self.f.pack()
        liste_exos_python = IntVar()
        liste_exos_python.set(1)  # ON
        Checkbutton(self.f, text="exos Python list", variable=liste_exos_python, command=self.liste_python).grid(row=14,
                                                                                                                 column=0)
        liste_exos_fortran = IntVar()
        liste_exos_fortran.set(0)  # ON
        Checkbutton(self.f, text="exos Fortran list", variable=liste_exos_fortran, command=self.liste_fortran).grid(
            row=15, column=0)

    def def_Langage(self):
        """
        to show an example of check box
        not used here
        """
        choice = ["Python", "Fortran"]
        var = StringVar()
        var.set("Python")
        Label(self.f, text="Langage choice: ", font=("Times", "11", "bold")).grid(row=17, column=0)
        OptionMenu(self.f, var, *choice).grid(row=7, column=1)

    def aide(self):
        """
        to display the help
        """

        def Quit():
            f_aide.destroy()

        """
        Help menu
        """
        f_aide = Tk()
        f_aide.geometry('700x200+40+40')
        F = Frame(f_aide)
        f_aide.title("Help")
        T = Text(f_aide, height=6, width=140, font=("Helvetica", "12", "italic"))
        T.pack()
        T.insert(END, self.aide_texte)
        Button(f_aide, text="Exit", bg="tomato", font=("Helvetica", "12", "bold"), command=Quit).pack(side="bottom")
        F.pack()
        f_aide.mainloop()

    def apropos(self):
        """ About box """
        txt = """ Fundamental aerodynamic suite : 
        Solutions found in the book : 
        'Exercices d'Aérodynamique Fondamentale',
        C. Airiau,
        A. Giovannini,
        P. Brancher
        Codes written by C. Airiau
        """
        self.root.option_add("*Dialog.msg.width", 35)
        showinfo("About", txt)

    def on_selected(self, x):
        """
        details of the selected item
        """
        self.selected_item = True
        self.s = self.tree.item(self.tree.focus())
        print("selection : ",self.s)

        # if s["tags"]!=["chap"] :
        #     print("codage  : ",s["tags"])
        #     print("commentaire : ",s["values"][2])
        #     if s["values"][1]=="f90":
        #         for k,liste in enumerate(exercise_list_fortran(display=False)):
        #             if s["values"][0] == liste['exo']:
        #                 print("valeur trouvée = ", k)

    def alert(self):
        """ Alert box """
        showinfo("alerte", "To implement")

    def browse_source(self):
        self.dirSrc.set(filedialog.askdirectory())

    def browse_run(self):
        self.dirRun.set(filedialog.askdirectory())

    def browse_runPy(self):
        self.dirRunPy.set(filedialog.askdirectory())

    def Sortir(self, event=""):
        print("Quit")
        try:
            os.remove(str(Path.home()) + "/ToDelete.txt")
        except:
            pass
        finally:
            sys.exit()

    def _liste_python_(self):
        """
        utiliser en cas de liste à cocher, inutile pour le moment
        """
        print("afficher la liste des programmes en python")

    def _liste_fortran_(self):
        """
        utiliser en cas de liste à cocher, inutile pour le moment
        """
        print("afficher la liste des programmes en fortran")
