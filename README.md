# Battle of Garbage Colelctors
This paper con- ducted a comprehensive comparison of the GC algorithms used in Java and Golang, two widely adopted programming languages with distinct approaches to memory management. By porting the GC implementation from Golang(v1.4) to the Java(OpenJDK@21), we performed an apples-to-apples comparison, evaluating the performance of both GC systems using a set of well-designed benchmarks. The objective is to provide insights into the trade-offs, strengths, and weaknesses of the GC implementations in these two languages, enabling developers to make informed decisions when choosing the appropriate language and the underlying Garbage Collection algorithm depending on their specific workload.

## Garbage Collection Algorithms
We chose SerialGC, ParallelGC, EspilonGC and Garbage First(G1)GC as the algorithms to compare our implementation of Go GC.

## Evaulation
We evaulated the garbage colelctors on 5 different benchmarks, including network-sentitive, memory-sentitive and IO-sentitive applications like Appache Kafka, H2 Database and GraphChi. Appache Kafka is a distributed event streaming platform which requires low latency for optimal performance, H2 is a memory intesive database and GraphiChi is a disk intesive workload. This makes these benchmarks get affected by Garbage Collection algorithms.

# GoGC implementation!

For build instructions please see the
[online documentation](https://openjdk.org/groups/build/doc/building.html),
or either of these files:

- [doc/building.html](doc/building.html) (html version)
- [doc/building.md](doc/building.md) (markdown version)

See <https://openjdk.org/> for more information about the OpenJDK
Community and the JDK and see <https://bugs.openjdk.org> for JDK issue
tracking.
