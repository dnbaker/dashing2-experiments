If you see errors about verifying signatures during building of the image, then clear your various Docker caches:

```
docker image purge
docker container purge
docker builder purge
```
