define(['jquery', 'app/ui/layout', 'bootstrap'], function($) {

    $(function() {

        // FAQ
        var $faq = $('#faq');
        if ($faq.length) {
            $faq.on('click', 'a', function(e) {
                var target = $(this).attr('href');
                $faq.find('.collapse').not(target).collapse('hide');
                $(target).collapse('show');

                e.preventDefault();
                return false;
            });
        }
    });
});
