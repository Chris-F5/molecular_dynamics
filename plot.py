import sys
import matplotlib.pyplot as plt
import events

frames = events.load_events(sys.stdin)
frames['energy']['total'] = frames['energy'][['kenetic', 'potential']].sum(axis=1)
print(frames['energy'])

fig, pos_ax = plt.subplots()
e_ax = pos_ax.twinx()
e_ax.axhline(y=0, color='black')
frames['pos'].loc[frames['pos']['id'] == 0].plot(y='x', kind='line', ax=pos_ax)
frames['pos'].loc[frames['pos']['id'] == 1].plot(y='x', kind='line', ax=pos_ax)
frames['energy'].plot(y='kenetic', kind='line', ax=e_ax, color='green')
frames['energy'].plot(y='potential', kind='line', ax=e_ax, color='purple')
frames['energy'].plot(y='total', kind='line', ax=e_ax, color='red')
plt.legend()
plt.show()
