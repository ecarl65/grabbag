import PySimpleGUI as sg

layout = [ 
        [sg.Text('Heading', key='-TEXT-')],
        [sg.Button('Read', key='-BUTTON-'), sg.Input('Read Input', tooltip="Make it a comma separated list", key='-READ-IN-')],
        [sg.Button('Write', key='-WRITE-'), sg.Input('Write Input', key='-WRITE-IN-')]
]

window = sg.Window('My new window', layout, finalize=True)

window['-TEXT-'].update('My new text value')

while True:             # Event Loop
    event, values = window.read()
    if event == sg.WIN_CLOSED:
        break
    elif event == "-BUTTON-":
        print(values['-READ-IN-'])
        window['-BUTTON-'].Update(disabled=True)
    elif event == "-WRITE-":
        print(values['-WRITE-IN-'])
        if window['-BUTTON-'].Disabled:
            print("Not printing read box, disabled")
        else:
            print(values['-READ-IN-'])

