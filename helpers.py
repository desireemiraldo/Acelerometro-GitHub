import Tkinter, tkFileDialog, helpers

#opens save file dialog
def file_save(text2save):
    f = tkFileDialog.asksaveasfile(mode='w', defaultextension=".txt")
    if f is None:
        return False
    f.writelines(text2save)
    f.close()
    return True

# opens dialog to open specific file
def file_open():
    file_path = tkFileDialog.askopenfilename()
    return file_path

