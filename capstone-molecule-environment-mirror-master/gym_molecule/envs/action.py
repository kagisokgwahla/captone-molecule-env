class Action():
    """
    This class is used to handle the actions taken in the environment
    """
    def __init__(self):
        """
        Creates the variables associated with the class
        """
        self.action_c = ''
        self.pos = ''   #front or back
        self.mol = ''
        self.query = ''
        self.isSmarts  = False

    def setAction(self,action,pos='front',query='',mol='C',isSmarts=False): #mol
        """
        Sets values for the variables associated with the class

        :type action: string
        :param action: The type of action

        :type pos: string, optional
        :param pos: The position where the action should be taken, defaults to front

        :type query: numpy.Array
        :param query: An array of molecule features to be removed from the molecule

        :type mol: string, optional
        :param mol: The molecule associated with the action to be taken, defaults to C
        """
        self.action_c  = action
        self.mol       = mol
        self.pos       = pos 
        self.query     = query
        self.isSmarts  = isSmarts       