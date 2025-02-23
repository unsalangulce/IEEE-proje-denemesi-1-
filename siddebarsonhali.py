
import tkinter as tk

root = tk.Tk()
root.geometry('400x600')
root.title('Single-Cell UMAP Viewer')
menu_bar_colour = '#383838'

toggle_icon = tk.PhotoImage(file='toggle.png')
home_icon = tk.PhotoImage(file='home.png')
service_icon = tk.PhotoImage(file='service.png')
update_icon = tk.PhotoImage(file='update.png')
contact_icon = tk.PhotoImage(file='contact.png')
about_icon = tk.PhotoImage(file='about.png')
close_icon = tk.PhotoImage(file='close.png')

def switch_indication(indicator_lb):
    home_btn_indicator,config(bg=menu_bar_colour)
    service_btn_indicator,config(bg=menu_bar_colour)
    contact_btn_indicator,config(bg=menu_bar_colour)
    update_btn_indicator,config(bg=menu_bar_colour)
    about_btn_indicator,config(bg=menu_bar_colour)
    
    indicator_lb.configure(bg='white')

menu_bar_frame = tk.Frame(root, bg=menu_bar_colour)
menu_bar_frame.pack(side=tk.LEFT, fill=tk.Y, pady=4, padx=3)
menu_bar_frame.pack_propagate(False)
menu_bar_frame.configure(width=50)

home_btn = tk.Button(menu_bar_frame, image=home_icon, bg=menu_bar_colour, bd=0, activebackground=menu_bar_colour,
                     command=lambda: switch_indication(indicatoe_lb=home_btn_indicator))
home_btn.place(x=9, y=70, width=30, height=40)

home_btn_indicator = tk.Label(menu_bar_frame, bg='white')
home_btn_indicator.place(x=3, y=70, width=5, height=40)

service_btn = tk.Button(menu_bar_frame, image=service_icon, bg=menu_bar_colour, bd=0, activebackground=menu_bar_colour)
service_btn.place(x=9, y=130, width=30, height=40)

service_btn_indicator = tk.Label(menu_bar_frame, bg='white')
service_btn_indicator.place(x=3, y=130, width=5, height=40)

update_btn = tk.Button(menu_bar_frame, image=update_icon, bg=menu_bar_colour, bd=0, activebackground=menu_bar_colour,
                       command=lambda: switch_indication(indicate_lb=update_btn_indicator))
update_btn.place(x=9, y=190, width=30, height=40)

update_btn_indicator = tk.Label(menu_bar_frame, bg='white')
update_btn_indicator.place(x=3, y=190, width=5, height=40)

contact_btn = tk.Button(menu_bar_frame, image=contact_icon, bg=menu_bar_colour, bd=0, activebackground=menu_bar_colour,
                        command=lambda: switch_indication(indicator_lb=contact_btn_indicator))
contact_btn.place(x=9, y=250, width=30, height=40)

contact_btn_indicator = tk.Label(menu_bar_frame, bg='white')
contact_btn_indicator.place(x=3, y=250, width=5, height=40)

about_btn = tk.Button(menu_bar_frame, image=about_icon, bg=menu_bar_colour, bd=0, activebackground=menu_bar_colour)
about_btn.place(x=9, y=310, width=30, height=40)

about_btn_indicator = tk.Label(menu_bar_frame, bg='white')
about_btn_indicator.place(x=3, y=310, width=5, height=40)

toggle_menu_btn = tk.Button(menu_bar_frame, image=toggle_icon, bg=menu_bar_colour, bd=0, activebackground=menu_bar_colour)
toggle_menu_btn.place(x=3, y=10, width=30, height=40)

root.mainloop()
