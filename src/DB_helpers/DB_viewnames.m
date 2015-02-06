function t = DB_viewnames(DBconn)
% t = DB_viewnames(DBconn)

t = fetch(DBconn, 'select name from sqlite_master where type is "view"');
