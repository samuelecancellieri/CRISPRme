import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html

import flask
import pandas as pd
import time
import os
from pages import navbar_creation

server = flask.Flask('app')

app = dash.Dash('app', server=server)

navbar = navbar_creation.Navbar()
app.layout = html.Div(
    [
        navbar,
        html.P('We are sorry, CRISPRme website is under maintenance, we will be back online as soon as possible.\nThanks for yout patience.')
    ]
)
if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=80, debug=False,
                   dev_tools_ui=False, dev_tools_props_check=False)
