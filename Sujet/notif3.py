#!/usr/bin/env python3
import os,subprocess,notify2,pyinotify
from gi.repository.GLib import GError

def auto_compile(path):
	wm = pyinotify.WatchManager()
	notifier = pyinotify.Notifier(wm)
	wm.add_watch(path, pyinotify.IN_CLOSE_WRITE,onSave)
	notifier.loop()
	
def onSave(event):
	if ".tex" in event.name:
		notify2.init("Mon appli python")
		try:
			subprocess.check_call(['mkdir', '-p', 'tmp'])
			subprocess.check_call(['latexmk', '-pdf', '-output-directory=tmp', '-shell-escape' ,'-halt-on-error' ,'Sujet_TP.tex'])
			subprocess.check_call(['mv', 'tmp/Sujet_TP.pdf', '.'])
			#subprocess.check_call(['mv', 'main.aux', 'tmp'])
			#subprocess.check_call(['mv', 'main.fdb_latexmk', 'tmp'])
			#subprocess.check_call(['mv', 'main.fls', 'tmp'])
			#subprocess.check_call(['mv', 'main.log', 'tmp'])
			#subprocess.check_call(['mv', 'main.nav', 'tmp'])
			#subprocess.check_call(['mv', 'main.out', 'tmp'])
			#subprocess.check_call(['mv', 'main.snm', 'tmp'])
			#subprocess.check_call(['mv', 'main.toc', 'tmp'])
		except subprocess.CalledProcessError:
			message = notify2.Notification("Build failed.")
		else:
			message = notify2.Notification("Build done.")
		try:
			message.show()
		except GError:
			print('Gerror')

def main():
	auto_compile([x[0] for x in os.walk(os.getcwd())])

if __name__ == '__main__':
    main()


