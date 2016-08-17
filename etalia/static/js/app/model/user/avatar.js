define(['app', 'app/model/user/user'], function (App) {

    App.Model.Avatar = App.Backbone.RelationalModel.extend({
        //urlRoot: App.config.api_root + '/user/users',

        defaults: {
            primary: null,
            avatar: null,
            date_uploaded: null,
            user: null
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'user',
                relatedModel: App.Model.User,
                includeInJSON: false,
                autoFetch: false
            }
        ]
    });

    /**
     * Handlebars helpers.
     */
    App.Handlebars.registerHelper('user_avatar', function(user) {
        if (!user) {
            return 'Expected user as first argument';
        }
        var avatar = null,
            avatars = App.getProperty(user, 'avatars');
        if (avatars.length) {
            avatar = avatars.find(function(avatar) {
                return avatar.get('primary');
            });
        }
        if (!avatar) {
            avatar = {
                avatar: '/static/img/avatar.jpg'
            }
        }
        return new App.Handlebars.SafeString(
            '<img src="' + App.getProperty(avatar, 'avatar') + '" alt="' + App.getProperty(user, 'first_name') + ' ' +
                App.getProperty(user, 'last_name') + '" width="200" height="200">'
        );
    });

    return App.Model.Avatar;
});
