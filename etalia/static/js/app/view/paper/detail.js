define([
    'app',
    'text!app/templates/paper/detail.hbs',
    'app/view/ui/modal',
    'app/view/paper/neighbors',
    'app/view/paper/related-threads',
    'altmetric'
], function (App, template) {

    App.View.Paper = App.View.Paper || {};

    var buttonsDefaults = {
        pin: false,
        ban: false,
        add: false,
        trash: false
    };

    return App.View.Paper.Detail = App.Backbone.View.extend({
        tagName: 'div',
        className: 'inner',

        template: App.Handlebars.compile(template),
        buttons: buttonsDefaults,

        events: {
            "click .detail-pin": "onPinClick",
            "click .detail-ban": "onBanClick",
            "click .detail-add": "onAddClick",
            "click .detail-trash": "onTrashClick",

            "click .detail-mail": "onMailClick",
            "click .detail-twitter": "onTwitterClick",
            "click .detail-google-plus": "onGooglePlusClick"
        },

        initialize: function (options) {
            if (options.buttons) {
                this.buttons = App._.extend(buttonsDefaults, options.buttons);
            }

            this.listenTo(this.model, "change", this.render);
            this.listenTo(this.model, "change:state", this.render);
        },

        onPinClick: function(e) {
            e.preventDefault();

            this.model.getState().togglePinned();
        },

        onBanClick: function(e) {
            e.preventDefault();

            this.model.getState().toggleBanned();
        },

        onAddClick: function(e) {
            e.preventDefault();

            this.model.getState().add();
        },

        onTrashClick: function(e) {
            e.preventDefault();

            this.model.getState().trash();
        },

        onMailClick: function(e) {
            e.preventDefault();

            window.location.href = 'mailto:'
                + '?subject=' + this.model.get('title')
                + '&body=Hi,I found this article and thought you might like it: '
                + this.model.get('url');

            App.trigger('etalia.paper.share', this.model, 'mail');
        },

        onTwitterClick: function(e) {
            e.preventDefault();

            var longURL = this.model.get('url'),
                title = this.model.get('title'); // TODO + main author
            App.$.ajax({
                type: 'POST',
                url: "https://www.googleapis.com/urlshortener/v1/url?key=AIzaSyCG85OFeMEgZMeHOI2dJB4VkuP-2HfGPPo",
                data: JSON.stringify({ longUrl:longURL}),
                contentType: 'application/json; charset=utf-8',
                success: function(data) {
                    App.popup(
                        'https://twitter.com/intent/tweet/'
                        + '?text=' + title
                        + '&url=' + encodeURI(data.id)
                        + '&via=etaliaio'
                        //+ '&hashtags=web,development';
                        , 'share-popup',
                        520, 377
                    );
                }
            });

            App.trigger('etalia.paper.share', this.model, 'twitter');
        },

        onGooglePlusClick: function(e) {
            e.preventDefault();

            var url = 'https://plus.google.com/share'
                    + '?url=' + this.model.get('url');

            App.popup(url, 'share-popup');

            App.trigger('etalia.paper.share', this.model, 'google-plus');
        },

        render: function() {
            App.log('PaperDetailView::render', this.model.get('id'));

            var that = this,
                attributes = App._.extend({}, this.model.attributes, {
                    pin_button: this.buttons.pin,
                    ban_button: this.buttons.ban,
                    add_button: this.buttons.add && !this.model.isAdded(),
                    trash_button: this.buttons.trash && this.model.isAdded()
                });

            this.$el.html(this.template(attributes));

            // Related Threads
            this.pushSubView(
                App.View.Paper.RelatedThreads.create({
                    paper_id: this.model.get('id'),
                    buttons: this.buttons,
                    return_callback: function() {
                        App.trigger('etalia.navigate', '/papers/' + that.model.get('slug') + '/');
                    }
                }, {
                    $target: this.$('[data-related-threads-placeholder]')
                })
            );

            // Neighbors
            this.pushSubView(
                App.View.Paper.Neighbors.create({
                    paper_id: this.model.get('id'),
                    buttons: this.buttons,
                    return_callback: function() {
                        App.trigger('etalia.navigate', '/papers/' + that.model.get('slug') + '/');
                    }
                }, {
                    $target: this.$('[data-neighbors-placeholder]')
                })
            );

            this.trigger('rendered');

            return this;
        },

        postRender: function() {
            _altmetric_embed_init();
        }
    });
});
