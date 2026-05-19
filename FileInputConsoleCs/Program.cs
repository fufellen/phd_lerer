using System;
using System.IO;
using System.Text;

namespace FileInputConsoleCs
{
    internal static class Program
    {
        private static int Main(string[] args)
        {
            Console.OutputEncoding = Encoding.UTF8;

            string filePath;
            if (args.Length > 0)
            {
                filePath = args[0].Trim();
                Console.WriteLine("Input file: " + filePath);
            }
            else
            {
                Console.Write("Enter input file path: ");
                filePath = (Console.ReadLine() ?? string.Empty).Trim();
            }

            if (filePath.Length == 0)
            {
                Console.Error.WriteLine("Error: empty file path.");
                return 2;
            }

            string fullPath = Path.GetFullPath(filePath);
            if (!File.Exists(fullPath))
            {
                Console.Error.WriteLine("Error: file was not found: " + fullPath);
                return 3;
            }

            FileInfo info = new FileInfo(fullPath);
            Console.WriteLine("Opened: " + info.FullName);
            Console.WriteLine("Size: " + info.Length + " bytes");
            Console.WriteLine("----- file content -----");

            int lineNumber = 1;
            using (StreamReader reader = new StreamReader(fullPath, Encoding.UTF8, true))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    Console.WriteLine("{0,4} | {1}", lineNumber, line);
                    lineNumber++;
                }
            }

            Console.WriteLine("----- end -----");
            Console.WriteLine("Lines: " + (lineNumber - 1));

            return 0;
        }
    }
}
