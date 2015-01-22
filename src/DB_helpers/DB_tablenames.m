function t = DB_tablenames(DBconn)

t = fetch(DBconn, 'select name from sqlite_master where type is "table"');
