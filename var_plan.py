# Make print work the same in all versions of Python, set up numpy,
# matplotlib, and use a nicer set of plot parameters:
# from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import datetime
import pytz
import sys
import matplotlib.dates as mdates
from astropy.coordinates import get_moon, get_sun
# from astropy.utils import conf


plt.style.use(astropy_mpl_style)
print sys.getrecursionlimit()
sys.setrecursionlimit(5500)


def read_cat(filename):
    f = open(filename)
    statr_arraay = []
    names, mags, Ps = [], [], []
    for line in f:
        if line[0] != "#":
            ll = line.split()
            RA = ll[0] + ":" + ll[1] + ":" + ll[2]
            DEC = ll[3] + ":" + ll[4] + ":" + ll[5]
            name = line[24:34]
            mag = line[72:78]
            P = line[83:89]
            # print name
            star = SkyCoord(RA, DEC, frame='icrs', unit=(u.hourangle, u.deg))
            statr_arraay.append(star)
            names.append(name)
            mags.append(mag)
            Ps.append(P)

    return statr_arraay, names, mags, Ps


def plot_star(i, location, frame_this_to_next_night, time_plot, time, sunaltazs_this_to_next_night, moonaltazs_this_to_next_night, slen):
    plt.clf()
    star = stars[i]
    star_name = star_names[i]

    # star_altaz = star.transform_to(AltAz(obstime=time, location=Derenivka))
    # print("Star's Altitude = {0.alt:.5} on 23:00:00".format(star_altaz))

    # star_altazs_this_night = star.transform_to(frame_this_night)

    star_altazs_this_to_next_night = star.transform_to(frame_this_to_next_night)

    ##############################################################################
    # Make a beautiful figure illustrating nighttime and the altitudes of M33 and
    # the Sun over that time:
    # print delta_midnight
    # fig, ax = plt.subplots()

    plt.title(str(i) + "/" + str(slen) + " " + star_name + ' [mag=' + mags[i] + ', P=' + str(float(Ps[i]) * 24.0) + ' h ]')
    plt.scatter(time_plot, star_altazs_this_to_next_night.alt.value, label=star_name,
                c=star_altazs_this_to_next_night.az.value, lw=0, s=8, cmap='viridis')

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))  # "%Y-%m-%d %H:%M"
    plt.fill_between(time_plot, 0, 90,
                     sunaltazs_this_to_next_night.alt < -0 * u.deg, color='0.5', zorder=0)
    plt.fill_between(time_plot, 0, 90,
                     sunaltazs_this_to_next_night.alt < -18 * u.deg, color='k', zorder=0)

    plt.plot(time_plot, sunaltazs_this_to_next_night.alt, color='r', label='Sun')
    plt.plot(time_plot, moonaltazs_this_to_next_night.alt, color='y', ls='--', label='Moon')  # color - [0.75]*3

    plt.colorbar().set_label('Azimuth [deg]')
    plt.legend(loc='upper left')
    # plt.xlim(-12, 12)
    plt.xlim(min(time_plot), max(time_plot))
    plt.ylim(0, 90)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Altitude [deg]')
    plt.gcf().autofmt_xdate()
    plt.show()


def on_press(event):
    global i
    # print('you pressed button No', event.button)
    if event.button == 1:
        i = i - 1
        if i < 0:
            i = 0
        plot_star(i, Derenivka, frame_this_to_next_night, time_plot, time, sunaltazs_this_to_next_night, moonaltazs_this_to_next_night, slen)
    elif event.button == 3:
        i = i + 1
        plot_star(i, Derenivka, frame_this_to_next_night, time_plot, time, sunaltazs_this_to_next_night, moonaltazs_this_to_next_night, slen)


stars, star_names, mags, Ps = read_cat("lkd-VAR_test.txt")
slen = len(stars)

##############################################################################
# Use `astropy.coordinates.EarthLocation` to provide the location of Bear
# Mountain and set the time to 11pm EDT on 2012 July 12:
Derenivka = EarthLocation(lat=48.63 * u.deg, lon=22.53 * u.deg, height=385 * u.m)

# utcoffset = +2*u.hour  # Eastern Daylight Time !!!!!!!!!!!!!!!!!!!!!!!!!!

now = datetime.datetime.now(pytz.timezone('Europe/Kiev'))
now_date = now.strftime("%Y-%m-%d")
utcoffset = now.utcoffset().total_seconds() / 60 / 60
utcoffset = utcoffset * u.hour

time = Time(now_date + ' 23:00:00') - utcoffset
# print time

# sys.exit()

#########################################
midnight = Time(now_date + ' 00:00:00') - utcoffset
delta_midnight = np.linspace(-2, 10, 100) * u.hour
# print delta_midnight
frame_this_night = AltAz(obstime=midnight + delta_midnight, location=Derenivka)

delta_midnight = np.linspace(-8, 10, 100) * u.hour
times_this_to_next_night = midnight + delta_midnight

time_plot = []
for t in times_this_to_next_night:

    t.format = 'datetime'
    tt = t.value
    time_plot.append(tt)


frame_this_to_next_night = AltAz(obstime=times_this_to_next_night, location=Derenivka)
sunaltazs_this_to_next_night = get_sun(times_this_to_next_night).transform_to(frame_this_to_next_night)
moon_this_to_next_night = get_moon(times_this_to_next_night)
moonaltazs_this_to_next_night = moon_this_to_next_night.transform_to(frame_this_to_next_night)
#########################################


global i
i = 0
fig, ax = plt.subplots()
cid = fig.canvas.mpl_connect('button_press_event', on_press)
# connection_id = fig.canvas.mpl_connect('button_press_event', onclick)
# fig.canvas.mpl_connect('pick_event', onpick)
# plt.tight_layout()

plot_star(i, Derenivka, frame_this_to_next_night, time_plot, time, sunaltazs_this_to_next_night, moonaltazs_this_to_next_night, slen)
