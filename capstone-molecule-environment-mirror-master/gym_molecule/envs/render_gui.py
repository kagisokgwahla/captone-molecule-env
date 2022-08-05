from tkinter import *
from PIL import ImageTk, Image


class Render():
    """
    This class handles rendering of the environment
    """
    
    def __init__(self, master):
        """
        This is the constructor

        :param master: The tkinter window to be used for the GUI
        :type master: tkinter.Toplevel
        """
        self.mols = []
        self.obs = []
        self.frames = []
        
        self.index = -1
        self.current = -1

        self.topFrame = Frame(master)
        self.topFrame.pack()
        self.bottomFrame = Frame(master)
        self.bottomFrame.pack(side = BOTTOM)
        self.midFrame = Frame(master)
        self.midFrame.pack(side = BOTTOM)

        self.molLabel = Label(self.topFrame)
        self.molLabel.pack()
        self.graph = Label(self.midFrame)
        self.graph.pack()

    def update(self, mol):
        """
        Updates information need for render

        :param mol: Image of the molecule
        :type mol: PIL.ImageTk
        """
        self.mols.append(mol)
        self.index += 1

        for i in range(1,360):
            pic = "gym_molecule/envs/resources/"+str(i)+".png"
            image = Image.open(pic).resize((400,400), Image.ANTIALIAS)
            frame = ImageTk.PhotoImage(image)
            self.frames.append(frame)
        self.obs.append(self.frames)

    def updateGif(self, ind):
        """
        changes the graph image in order to animate the graph's rotation

        :param ind: The index of the image to be displayed
        :type ind: int
        """
        frame = self.frames[ind]
        if ind == 358:
            ind = 0
        else:
            ind +=1
        self.graph.configure(image = frame)
        self.graph.image = frame
        self.midFrame.after(40, self.updateGif, ind)

    def nxt(self):
        """
        Iterates to the next step
        """
        if self.current < self.index:
            self.current += 1
            pic = ImageTk.PhotoImage(self.mols[self.current])
            self.molLabel.configure(image=pic)
            self.molLabel.image = pic

            self.frames = self.obs[self.current]
            self.midFrame.after(0, self.updateGif, 0)

    def prev(self):
        """
        Iterates to the previous step
        """
        if self.current > 0:
            self.current -= 1
            pic = ImageTk.PhotoImage(self.mols[self.current], master = self.molLabel)
            self.molLabel.configure(image=pic)
            self.molLabel.image = pic

            self.frames = self.obs[self.current]
            self.midFrame.after(0, self.updateGif, 0)

    def render(self):
        """
        Diplays the state
        """

        if len(self.mols) != 0:
            self.current = self.index
            pic = ImageTk.PhotoImage(self.mols[self.index])
            self.molLabel.configure(image=pic)
            self.molLabel.image = pic

            prevButton = Button(self.bottomFrame, text="PREVIOUS", command = self.prev)
            prevButton.pack(side=LEFT)
            nextButton = Button(self.bottomFrame, text="NEXT", command = self.nxt)
            nextButton.pack(side=RIGHT)
            self.midFrame.after(0, self.updateGif, 0)

    def reset(self):
        """
        resets the render information
        """
        self.mols = []
        self.obs = []
        self.frames = []

