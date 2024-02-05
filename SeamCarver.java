import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.StdOut;

public class SeamCarver {
    private static Picture picCopy;
    private int width;
    private int height;
    private double[][] distTo;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        picCopy = new Picture(picture);
        width = picture.width();
        height = picture.height();

        distTo = new double[width][height];
    }

    // current picture
    public Picture picture() {
        return new Picture(picCopy);
    }

    // width of current picture
    public int width() {
        return width;
    }

    // height of current picture
    public int height() {
        return height;
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || x > width - 1 || y < 0 || y > height - 1)
            throw new IllegalArgumentException("Invalid coordinates");
        int up = y - 1;
        int down = y + 1;
        int left = x - 1;
        int right = x + 1;
        if (x == 0) {
            left = width - 1;
        }
        if (x == width - 1) {
            right = 0;
        }
        if (y == 0) {
            up = height - 1;
        }
        if (y == height - 1) {
            down = 0;
        }

        right = picCopy.getRGB(right, y);
        left = picCopy.getRGB(left, y);
        up = picCopy.getRGB(x, up);
        down = picCopy.getRGB(x, down);

        return Math.sqrt(gradient(left, right) + gradient(up, down));
    }

    // helper to calculate gradient
    private double gradient(int rgb1, int rgb2) {
        // calculated using the picture documentation
        int r1 = (rgb1 >> 16) & 0xFF;
        int r2 = (rgb2 >> 16) & 0xFF;
        int g1 = (rgb1 >> 8) & 0xFF;
        int g2 = (rgb2 >> 8) & 0xFF;
        int b1 = (rgb1) & 0xFF;
        int b2 = (rgb2) & 0xFF;

        return Math.pow(r1 - r2, 2) + Math.pow(g1 - g2, 2)
                + Math.pow(b1 - b2, 2);
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        // transpose picture
        transpose(picCopy);
        // StdOut.println("width: " + width);
        // StdOut.println("height: " + height);

        // save output
        int[] output = findVerticalSeam();

        // transpose it back
        transpose(picCopy);
        // StdOut.println("width: " + width);
        // StdOut.println("height: " + height);
        return output;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        double[][] energy = new double[width][height];
        // double[][] distTo = new double[width][height];

        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                energy[c][r] = energy(c, r);
            }
        }

        // initialize top row
        for (int c = 0; c < width; c++) {
            distTo[c][0] = energy(c, 0);
        }

        int[][] edgeTo = new int[width][height];

        // setting up edgeTo and distTo points
        for (int r = 1; r < height; r++) {
            for (int c = 0; c < width; c++) {
                // straight up
                double minChamp = distTo[c][r - 1] + energy[c][r];
                int idxChamp = c;
                if (c != 0) {
                    double topLCheck = distTo[c - 1][r - 1] + energy[c][r];
                    if (topLCheck < minChamp) {
                        minChamp = topLCheck;
                        idxChamp = c - 1;
                    }
                }
                if (c != width - 1) {
                    double topRCheck = distTo[c + 1][r - 1] + energy[c][r];
                    if (topRCheck < minChamp) {
                        minChamp = topRCheck;
                        idxChamp = c + 1;
                    }
                }
                distTo[c][r] = minChamp;
                edgeTo[c][r] = idxChamp;
            }
        }

        int idxChamp = 0;
        double distToChamp = Double.POSITIVE_INFINITY;
        for (int i = 0; i < width; i++) {
            double dist = distTo[i][height - 1];
            if (distToChamp > dist) {
                distToChamp = dist;
                idxChamp = i;
            }
        }

        int[] seam = new int[height];
        seam[height - 1] = idxChamp;
        for (int r = height - 2; r >= 0; r--) {
            seam[r] = edgeTo[seam[r + 1]][r + 1];
        }

        return seam;
    }

    // transpose picture
    private void transpose(Picture pic) {
        width = pic.height(); // set orig height to new width
        height = pic.width(); // set orig width to new height

        picCopy = new Picture(width, height);
        // based on orig set-up
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                // flip get and set coordinates
                picCopy.setRGB(c, r, pic.getRGB(r, c));
            }
        }
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (height == 1 || seam.length != width) {
            throw new IllegalArgumentException("wrong height ");
        }
        for (int i = 1; i < seam.length; i++) {
            int diff = seam[i] - seam[i - 1];
            if (!(diff >= -1 && diff <= 1)) {
                throw new IllegalArgumentException();
            }
        }
        transpose(picCopy);
        removeVerticalSeam(seam);
        transpose(picCopy);
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (width == 1 || seam.length != height) {
            throw new IllegalArgumentException();
        }
        for (int i = 1; i < seam.length; i++) {
            int diff = seam[i] - seam[i - 1];
            if (!(diff >= -1 && diff <= 1)) {
                throw new IllegalArgumentException();
            }
        }

        Picture newPicCopy = new Picture(width - 1, height);

        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width - 1; c++) {
                if (seam[r] < 0 || seam[r] > width - 1) {
                    throw new IllegalArgumentException();
                }
                if (c >= seam[r]) {
                    newPicCopy.setRGB(c, r, picCopy.getRGB(c + 1, r));
                }
                else {
                    newPicCopy.setRGB(c, r, picCopy.getRGB(c, r));
                }
            }
        }
        width--;
        picCopy = newPicCopy;
    }

    //  unit testing (required)
    public static void main(String[] args) {
        Picture picture = new Picture(args[0]);

        SeamCarver sc = new SeamCarver(picture);
        sc.findVerticalSeam();

        for (int r = 0; r < sc.height; r++) {
            for (int c = 0; c < sc.width; c++) {
                StdOut.printf("%7.2f ", sc.distTo[c][r]);
            }
            StdOut.println();
        }
    }

}
