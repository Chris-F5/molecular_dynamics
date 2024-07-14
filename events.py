import json
import pandas as pd

event_types = {
  'pos': ['time', 'id', 'x', 'y'],
  'force': ['time', 'id', 'x', 'y'],
  'energy': ['time', 'kenetic', 'potential'],
}

def load_events(f):
  events_by_type = {}
  for line in f:
    event = json.loads(line)
    type = event['type']
    del event['type']
    events_by_type[type] = events_by_type.get(type, []) + [event]

  frames = {type:pd.DataFrame(events) for type,events in events_by_type.items()}
  for type,frame in frames.items():
    frame.set_index(event_types[type][0], inplace=True)
  return frames
