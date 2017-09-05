from helpers import *

def main():
    print('Choose your csv data file to read from...')

    #Get path to target file
    filepath = file_open()
    index = 0
    linesNumber  = 31512
    
    #Opens target file
    with open(filepath) as file:

        #Read lines in file
        lines = file.readlines()
        while True:

             #Get next 31512 lines
            choosenLines = lines[((linesNumber * index) + 1) : (linesNumber * (index + 1))]

            #Checks if linex exist in file, if not program terminates
            if not choosenLines:
                return

            #Save file using dialog
            if not file_save(choosenLines):
                return

            index = index + 1
        print('Closing file...')

if __name__ == "__main__":
    main()
    




