#
# Regular cron jobs for the libspecbleach package.
#
0 4	* * *	root	[ -x /usr/bin/libspecbleach_maintenance ] && /usr/bin/libspecbleach_maintenance
