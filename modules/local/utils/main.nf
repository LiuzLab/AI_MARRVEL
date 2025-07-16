def showVersion() {
    if (!params.version) {
        return
    }

    println workflow.manifest.version
    exit 0
}
