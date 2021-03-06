# ParallelProgrammingCourse

`В процессе выполнения заданий использовались Linux Red Hat, Debian GNU/Linux`

## Практическое задание № 1 
### Задание
**Реализовать последовательный алгоритм матричного умножения и оценить влияние кэша на время выполнение программы(в зависимости от перестановки циклов).**  

Реализация задания и проведение анализа связанного с кэш памятью проводились на ВК **IBM POLUS**.  
Для компиляциии всех необходимых исходников и запуска исполняемых файлов подготовлен **Makefile**.

#### Описание .h и .cpp файлов:
 * **Matrix.h** - заголовочный файл, который описывает матрицу
 * **genMatrix.cpp** - генерация рандомых матриц в бинарном виде
 * **printMtr.cpp** - вывод содержимого файла с матрицами на экран
 * **createC0.cpp** - независимое вычисление произведения матриц с помощью библиотеки "ublas"
 * **main.cpp** - функция умножения матриц в следующих режимах: 0 - *ijk*, 1 - *ikj*, 2 - *kij*, 4 - *jki*, 5 - *kji* 
 * **compare.cpp** - сравнение полученного произведения матриц в результате умножения с заданным режимом и независимого умножения (ublas)

В результате практического задания формируется отчет, содержащий описание проделанной работы, графики работы программы (с помощью .sh - скриптов: **gnuScript.sh** и **plotScript.sh**)

## Практическое задание № 2
### Задание
**Реализовать последовательный алгоритм блочного матричного умножения и оценка влияния кэша на время выполнение программы(в зависимости от перестановки циклов).**  

Для выполнения этого практического задания используется (PAPI)
**Performance Application Programming Interface** — переносимый интерфейс, реализованный в виде библиотеки, для доступа к счетчикам аппаратной производительности на различных современных микропроцессорах ([wiki](https://ru.wikipedia.org/wiki/PAPI)).   
Реализация задания и проведение анализа связанного с кэш-памятью проводились на ВК **IBM POLUS**.  
Для компиляциии всех необходимых исходников и запуска исполняемых файлов подготовлен **Makefile**. 

#### Описание .h и .cpp файлов
 * **Matrix.h** - заголовочный файл, который описывает матрицу
 * **genMatrix.cpp** - генерация рандомых матриц в бинарном виде
 * **printMtr.cpp** - вывод содержимого файла с матрицами на экран
 * **createC0.cpp** - независимое вычисление произведения матриц с помощью библиотеки "ublas"
 * **PPT_prac2.cpp** - функция умножения матриц в следующих режимах: 0 - *ijk* (размер блока - `32x32`), 1 - *ikj* (размер блока - `32x32`),2 - *ikj* (`оптимальный` размер блока ```int blocksize = sqrt(cachesize/(3*sizeof(float)));```)
 * **compare.cpp** - сравнение полученного произведения матриц в результате умножения с заданным режимом и независимого умножения (ublas)

В результате практического задания формируется отчет, содержащий описание проделанной работы, графики работы программы (с помощью .sh - скриптов: **gnuScript.sh** и **plotScript.sh**

## Практическое задание № 3 
### Задание
__Реализация параллельного алгоритма поиска простых чисел в заданном диапазоне с помощью "решета Эратосфена".__

Разработка и тестирование производились на ВК **IBM Blue Gene/P**.

Для компиляции - `mpicxx main.cpp -o main`.

Примером для запуска на **IBM Blue GENE/P** - `mpisubmit.bg -n 3 -w 00:05:00 -m dual main -- 1 100 test_out.txt`.

В рамках этого задания была выполнена реализация алгоритма поиска простых чисел, с помощью `POSIX Threads`(PracticalTask3_5).

## Практическое задание № 4 
### Задание
__Реализовать параллельный алгоритм блочного умножения матриц, предусмотрев равномерное распределение матриц блоками строк по процессам или блоками столбцов.__  
*Дополнительным условием явялется исследование влияния мэппинга*.  
Реализация задания и проведение анализа связанного с кэш-памятью проводились на ВК **IBM Blue Gene/P**.  
Для компиляции всех необходимых исходников и запуска исполняемых файлов подготовлен **Makefile**.

#### Описание .h и .cpp файлов:
 * **Matrix.h** - заголовочный файл, который описывает матрицу
 * **genMatrix.cpp** - генерация рандомых матриц в бинарном виде
 * **printMtr.cpp** - вывод содержимого файла с матрицами на экран
 * **main.cpp** - функция умножения матриц (в зависимости от размеров матриц возможно строковое распределение или распределение столбцами)
 * **genMap.cpp** - генерация мэп файла (организация произвольного мэппинга)

В результате практического задания формируется отчет, содержащий описание проделанной работы, графики работы программы (с помощью .sh - скрипта: **script.sh**.

## Практическое задание № 5 
### Задание
__Реализовать параллельный алгоритм DNS умножения матриц.__  
*Дополнительным условием явялется исследование влияния мэппинга*.
Реализация задания и проведение анализа связанного с кэш-памятью проводились на ВК **IBM Blue Gene/P**.  
Для компиляции всех необходимых исходников и запуска исполняемых файлов подготовлен **Makefile**.

#### Описание .h и .cpp файлов:
 * **Matrix.h** - заголовочный файл, который описывает матрицу
 * **genMatrix.cpp** - генерация рандомых матриц в бинарном виде
 * **printMtr.cpp** - вывод содержимого файла с матрицами на экран
 * **main.cpp** - функция умножения матриц, спользуя декартову решетку процессов, с предварительным распределением блоков матрицы
 * **genMap.cpp** - генерация мэп файла (организация произвольного мэппинга)

В результате практического задания формируется отчет, содержащие описание проделанной работы, графики работы программы (с помощью .sh - скрипта: **script.sh**.  

## Практическое задание
### Задание
**Поиск численного решения параболического уравнения.**  
В результате полученная программа должна была выполняться корректно на вычислительных системах  
(**IBM Blue Gene/P** и **Polus**) МГУ.  
Для компиляции всех необходимых исходников и запуска исполняемых файлов подготовлен **Makefile**.  

В качестве входных параметров выступают размер сетки, а также параметры дискретизации.
