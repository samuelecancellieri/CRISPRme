import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html

import flask
import pandas as pd
import time
import os

server = flask.Flask('app')

app = dash.Dash('app', server=server)

app.scripts.config.serve_locally = False
dcc._js_dist[0]['external_url'] = 'https://cdn.plot.ly/plotly-basic-latest.min.js'

app.layout = html.Div(
    'We are sorry, CRISPRme website is under maintenance, we will be back online as soon as possible.\nThanks for yout patience.')

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=80, debug=False,
                   dev_tools_ui=False, dev_tools_props_check=False)
