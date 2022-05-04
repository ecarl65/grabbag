import PySimpleGUIWeb as sg

layout = [  [sg.Text('Row 1'), sg.Text("What's your name?")],
            [sg.Text('Row 2'), sg.Input()],
            [sg.Text('Row 3'), sg.Button('Ok')] ]

event, values = sg.Window('List Comprehensions', layout).read(close=True)
