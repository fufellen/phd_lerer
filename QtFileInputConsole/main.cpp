#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QStringList>
#include <QTextStream>
#include <QtGlobal>

#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
#include <QStringConverter>
#else
#include <QTextCodec>
#endif

static void setupUtf8(QTextStream &stream)
{
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    stream.setEncoding(QStringConverter::Utf8);
#else
    stream.setCodec(QTextCodec::codecForName("UTF-8"));
#endif
}

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);
    app.setApplicationName("QtFileInputConsole");

    QTextStream input(stdin);
    QTextStream output(stdout);
    QTextStream error(stderr);
    setupUtf8(input);
    setupUtf8(output);
    setupUtf8(error);

    QString filePath;
    const QStringList args = app.arguments();

    if (args.size() > 1) {
        filePath = args.at(1).trimmed();
        output << "Input file: " << filePath << '\n';
    } else {
        output << "Enter input file path: ";
        output.flush();
        filePath = input.readLine().trimmed();
    }

    if (filePath.isEmpty()) {
        error << "Error: empty file path.\n";
        return 2;
    }

    QFileInfo fileInfo(filePath);
    if (fileInfo.isRelative()) {
        fileInfo.setFile(QDir::current().absoluteFilePath(filePath));
    }

    if (!fileInfo.exists() || !fileInfo.isFile()) {
        error << "Error: file was not found: "
              << QDir::toNativeSeparators(fileInfo.absoluteFilePath()) << '\n';
        return 3;
    }

    QFile file(fileInfo.absoluteFilePath());
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        error << "Error: cannot open file: " << file.errorString() << '\n';
        return 4;
    }

    output << "Opened: " << QDir::toNativeSeparators(fileInfo.absoluteFilePath()) << '\n';
    output << "Size: " << fileInfo.size() << " bytes\n";
    output << "----- file content -----\n";

    QTextStream fileStream(&file);
    setupUtf8(fileStream);

    int lineNumber = 1;
    while (!fileStream.atEnd()) {
        output << QString("%1 | %2")
                      .arg(lineNumber, 4)
                      .arg(fileStream.readLine())
               << '\n';
        ++lineNumber;
    }

    output << "----- end -----\n";
    output << "Lines: " << (lineNumber - 1) << '\n';

    return 0;
}
